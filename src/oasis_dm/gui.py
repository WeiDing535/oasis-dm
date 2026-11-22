"""Desktop GUI for the oasis_dm density-modification engine.

A native window application built on PySide6 (Qt 6).  Launch with::

    oasis-dm-gui             # after: pip install ".[gui]"
    python -m oasis_dm.gui   # equivalent

The engine always runs as a *subprocess* of ``python -m oasis_dm`` (via
QProcess, so it can be stopped cleanly and its log streams into the
window), and the command line at the bottom shows the exact CLI
invocation the Run button will perform — the GUI doubles as a gentle
introduction to the command line interface.

Why PySide6 and not Tkinter?  Tkinter ships with CPython in theory, but in
practice the Tcl/Tk bundled with system Pythons is frequently broken — on
this repository's reference machine (macOS 15, Xcode toolchain Python
3.9.6) calling ``tk.Tk()`` aborts the whole process.  PySide6 is a single
``pip install`` away and its macOS wheels were verified to work there.

The module imports fine without PySide6 installed (the Qt imports happen
inside :func:`build_app`), so the pure-logic helpers below stay testable
everywhere.
"""
from __future__ import annotations

import os
import re
import shlex
import sys
import time
from pathlib import Path

# --- pure logic (no Qt; unit-tested in tests/test_gui.py) ------------------


def read_mtz_summary(path):
    """Column labels + header facts for the column pickers, via gemmi."""
    import gemmi

    m = gemmi.read_mtz_file(str(path))
    d = m.make_d_array()
    finite = [float(x) for x in d if x == x]    # drop NaN, cast float32
    return {
        "labels": list(m.column_labels()),
        "spacegroup": m.spacegroup.hm,
        "nref": int(len(m.array)),
        "dmin": min(finite) if finite else float("nan"),
    }


def pick_default(labels, wanted, fallback=None):
    """First label in `wanted` present in the file, else `fallback`.

    ``fallback=""`` means "leave the column empty" and is honoured; when no
    fallback is given the first label in the file is used.
    """
    have = set(labels)
    for cand in wanted:
        if cand in have:
            return cand
    if fallback is not None:
        return fallback
    return labels[0] if labels else ""


def build_argv(st):
    """GUI state dict -> engine argv.

    Only options that differ from the engine defaults are emitted, so the
    previewed command line stays short and honest.
    """
    argv = ["--mtzin", st["mtzin"], "--mtzout", st["mtzout"]]
    if st.get("stats"):
        argv += ["--stats", st["stats"], "--tag", st.get("tag") or "dm"]
    # column mapping: only when it differs from the engine default
    for key, flag, default in (
        ("fo", "--fo", "FOBS"), ("sigf", "--sigf", "SIGFOBS"),
        ("phib", "--phib", "PHIB"), ("fom", "--fom", "FOM"),
        ("free", "--free", "FreeR"), ("truth", "--truth", "PHIC"),
    ):
        if st.get(key) and st[key] != default:
            argv += [flag, st[key]]
    if st.get("solc") is not None and abs(st["solc"] - 0.5) > 1e-9:
        argv += ["--solc", "%.3f" % st["solc"]]
    if st.get("cycles") is not None and st["cycles"] != 8:
        argv += ["--cycles", str(st["cycles"])]
    if not st.get("auto_stop", True):
        argv += ["--no-auto-stop"]
    if st.get("mask", "variance") != "variance":
        argv += ["--mask", st["mask"]]
    if st.get("hist_mode", "cell") != "cell":
        argv += ["--hist-mode", st["hist_mode"]]
    if st.get("free_metric", "cc") != "cc":
        argv += ["--free-metric", st["free_metric"]]
    if not st.get("extend", True):
        argv += ["--no-extend"]
    if st.get("aniso"):
        argv += ["--aniso"]
    if st.get("per_shell_d"):
        argv += ["--per-shell-d"]
    # NCS: exactly one of none / operator file / auto-detect
    ncs = st.get("ncs", "none")
    if ncs == "file":
        argv += ["--ncs", st.get("ncs_file", "") or ""]
    elif ncs == "auto":
        argv += ["--ncs-auto"]
    # histogram library
    if st.get("no_hist"):
        argv += ["--no-hist"]
    elif st.get("hist_lib"):
        argv += ["--hist-lib", st["hist_lib"]]
    if st.get("sigmaq_source", "blend") != "blend":
        argv += ["--sigmaq-source", st["sigmaq_source"]]
    return argv


def validate_state(st):
    """Human-readable problems with a run request, [] when fine."""
    problems = []
    if not st.get("mtzin") or not Path(st["mtzin"]).is_file():
        problems.append("choose an existing input MTZ file")
    if not st.get("mtzout"):
        problems.append("choose where to write the output MTZ")
    solc = st.get("solc")
    if solc is None or not 0.05 < solc < 0.95:
        problems.append("solvent content must be between 0.05 and 0.95 "
                        "(typically 0.35-0.65)")
    if st.get("ncs") == "file" and not st.get("ncs_file"):
        problems.append("choose the NCS operator file, or pick another NCS "
                        "mode")
    return problems


def command_string(st):
    return "oasis-dm " + " ".join(shlex.quote(a) for a in build_argv(st))


# engine per-cycle table line:  cyc  meanFOM  freeCC  freeR  ...
CYCLE_RE = re.compile(
    r"^\s*(\d+)\s+(-?[0-9.]+)\s+(-?[0-9.]+)\s+(-?[0-9.]+)")


def parse_cycle_line(line):
    """'( 1  0.71  0.65 ...)' -> [cycle, meanFOM, freeCC] or None."""
    m = CYCLE_RE.match(line)
    if not m:
        return None
    cyc, fom, cc = int(m.group(1)), float(m.group(2)), float(m.group(3))
    if cyc > 999 or not (0.0 <= fom <= 1.0) or abs(cc) > 1.0:
        return None      # not the statistics table (NaN columns etc.)
    return [cyc, fom, cc]


def _pythonpath_env():
    """Process env that finds oasis_dm from a source checkout."""
    env = dict(os.environ)
    here = Path(__file__).resolve().parent          # .../src/oasis_dm
    src = here.parent                                # .../src
    if (src / "oasis_dm" / "__init__.py").exists() and str(src) \
            not in env.get("PYTHONPATH", ""):
        env["PYTHONPATH"] = str(src) + os.pathsep + env.get("PYTHONPATH", "")
    return env


def qt_available():
    try:
        import PySide6  # noqa: F401
    except ImportError:
        return False
    return True


# --- the application --------------------------------------------------------
#
# Everything Qt lives inside build_window() / build_app() so importing this
# module stays side-effect-free on machines without PySide6 (CI).

COLUMNS = (
    ("fo", "Fo", ("FOBS", "FP", "F")),
    ("sigf", "σ(Fo)", ("SIGFOBS", "SIGFP", "SIGF")),
    ("phib", "PHIB", ("PHIB", "PHI")),
    ("fom", "FOM", ("FOM", "WT")),
    ("free", "free set", ("FreeR", "FREE")),
    ("truth", "truth", ("PHIC",)),
)


def build_app():
    """Create the QApplication and the main window (returns both)."""
    from PySide6.QtWidgets import QApplication

    app = QApplication.instance() or QApplication(sys.argv or ["oasis-dm"])
    win = build_window()
    return app, win


def build_window():
    from PySide6.QtCore import Qt
    from PySide6.QtGui import QColor, QFont, QPainter, QPen
    from PySide6.QtWidgets import (
        QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog, QFormLayout,
        QGroupBox, QHBoxLayout, QLabel, QLineEdit, QMainWindow,
        QPlainTextEdit, QPushButton, QRadioButton, QScrollArea, QSpinBox,
        QVBoxLayout, QWidget)

    from . import __version__

    class CcChart(QWidget):
        """free-set cc per cycle, drawn with QPainter (no chart library)."""

        def __init__(self):
            super().__init__()
            self.pts = []

        def set_cycles(self, cycles):
            self.pts = [(c[0], c[2]) for c in cycles]
            self.update()

        def paintEvent(self, _event):
            p = QPainter(self)
            p.setRenderHint(QPainter.Antialiasing)
            w, h = self.width(), self.height()
            p.fillRect(0, 0, w, h, QColor("white"))
            pad, pt_, pb, pr = 44, 12, 22, 14
            if not self.pts:
                p.setPen(QColor("#999999"))
                p.drawText(self.rect(), Qt.AlignCenter,
                           "free-set cc per cycle")
                return
            ys = [y for _x, y in self.pts] + [0.0]
            lo, hi = min(ys), max(ys)
            if hi - lo < 1e-6:
                hi = lo + 1e-6

            def xy(i, y):
                x = pad + (w - pad - pr) * (i / max(len(self.pts) - 1, 1))
                yy = h - pb - (h - pt_ - pb) * ((y - lo) / (hi - lo))
                return x, yy

            p.setPen(QPen(QColor("#bbbbbb"), 1))
            p.drawLine(pad - 6, h - pb, w - pr, h - pb)
            p.drawLine(pad - 6, h - pb, pad - 6, pt_)
            p.setPen(QColor("#888888"))
            p.drawText(2, pt_ + 14, "%.2f" % hi)
            p.drawText(2, h - pb - 4, "%.2f" % lo)
            p.setPen(QPen(QColor("#0088cc"), 2))
            for i in range(len(self.pts) - 1):
                x1, y1 = xy(i, self.pts[i][1])
                x2, y2 = xy(i + 1, self.pts[i + 1][1])
                p.drawLine(int(x1), int(y1), int(x2), int(y2))
            small = QFont()
            small.setPointSizeF(8)
            p.setFont(small)
            for i, (cyc, y) in enumerate(self.pts):
                x, yy = xy(i, y)
                p.setBrush(QColor("#0088cc"))
                p.setPen(Qt.NoPen)
                p.drawEllipse(int(x) - 3, int(yy) - 3, 6, 6)
                p.setPen(QColor("#888888"))
                p.drawText(int(x) - 4, h - 6, str(cyc))

    class Win(QMainWindow):
        def __init__(self):
            super().__init__()
            self.setWindowTitle("oasis-dm %s — density modification"
                                % __version__)
            self.resize(960, 880)
            self.cycles = []
            self.proc = None
            self._t0 = 0.0
            self._build_ui()

        # -- layout ------------------------------------------------------

        def _build_ui(self):
            central = QWidget()
            scroll = QScrollArea()
            scroll.setWidgetResizable(True)
            scroll.setWidget(central)
            self.setCentralWidget(scroll)
            root = QVBoxLayout(central)
            root.setContentsMargins(14, 14, 14, 14)
            root.setSpacing(10)

            # files
            g = QGroupBox("Files")
            f = QFormLayout(g)
            f.setFieldGrowthPolicy(QFormLayout.AllNonFixedFieldsGrow)
            self.e_in = QLineEdit()
            self.e_in.setMinimumWidth(420)
            self.e_in.textChanged.connect(
                lambda _t: self._debounce_input())
            b = QPushButton("Browse…")
            b.clicked.connect(self.browse_in)
            row = QHBoxLayout()
            row.addWidget(self.e_in, 1)
            row.addWidget(b)
            f.addRow("input MTZ", row)
            self.e_out = QLineEdit()
            b = QPushButton("Browse…")
            b.clicked.connect(self.browse_out)
            row = QHBoxLayout()
            row.addWidget(self.e_out, 1)
            row.addWidget(b)
            f.addRow("output MTZ", row)
            self.l_info = QLabel("no input selected yet")
            self.l_info.setStyleSheet("color: #64748b;")
            f.addRow("", self.l_info)
            root.addWidget(g)

            # columns
            g = QGroupBox("Columns (auto-detected from the input)")
            h = QHBoxLayout(g)
            self.cols = {}
            for key, label, _want in COLUMNS:
                box = QVBoxLayout()
                box.addWidget(QLabel(label))
                c = QComboBox()
                self.cols[key] = c
                box.addWidget(c)
                h.addLayout(box)
            root.addWidget(g)

            # parameters
            g = QGroupBox("Parameters")
            h = QHBoxLayout(g)
            box = QVBoxLayout()
            box.addWidget(QLabel("solvent content"))
            self.s_solc = QDoubleSpinBox()
            self.s_solc.setRange(0.05, 0.95)
            self.s_solc.setSingleStep(0.01)
            self.s_solc.setValue(0.5)
            box.addWidget(self.s_solc)
            h.addLayout(box)
            box = QVBoxLayout()
            box.addWidget(QLabel("max cycles"))
            self.s_cycles = QSpinBox()
            self.s_cycles.setRange(1, 100)
            self.s_cycles.setValue(8)
            box.addWidget(self.s_cycles)
            h.addLayout(box)
            self.cb_stop = QCheckBox("free-set early stop")
            self.cb_stop.setChecked(True)
            h.addWidget(self.cb_stop)
            h.addStretch(1)
            for attr, label, values in (
                    ("c_mask", "mask", ("variance", "wang")),
                    ("c_hist", "histogram target", ("cell", "protein")),
                    ("c_metric", "free metric", ("cc", "rfree"))):
                box = QVBoxLayout()
                box.addWidget(QLabel(label))
                c = QComboBox()
                c.addItems(list(values))
                setattr(self, attr, c)
                box.addWidget(c)
                h.addLayout(box)
            root.addWidget(g)

            # options
            g = QGroupBox("Options")
            h = QHBoxLayout(g)
            self.cb_extend = QCheckBox("phase extension")
            self.cb_extend.setChecked(True)
            self.cb_aniso = QCheckBox("anisotropic correction")
            self.cb_pershell = QCheckBox("per-shell D in output")
            h.addWidget(self.cb_extend)
            h.addWidget(self.cb_aniso)
            h.addWidget(self.cb_pershell)
            h.addStretch(1)
            root.addWidget(g)

            # NCS
            g = QGroupBox(
                "Non-crystallographic symmetry — expert option: averaging "
                "only pays with tight operators and several copies")
            v = QVBoxLayout(g)
            h = QHBoxLayout()
            self.r_ncs = {}
            for key, label in (("none", "none"), ("file", "operator file"),
                               ("auto", "auto-detect (free-set gated)")):
                r = QRadioButton(label)
                self.r_ncs[key] = r
                if key == "none":
                    r.setChecked(True)
                r.toggled.connect(self._sync_ncs)
                h.addWidget(r)
            h.addStretch(1)
            v.addLayout(h)
            f = QFormLayout()
            self.e_ncs = QLineEdit()
            self.e_ncs.setEnabled(False)
            b = QPushButton("Browse…")
            b.clicked.connect(self.browse_ncs)
            row = QHBoxLayout()
            row.addWidget(self.e_ncs, 1)
            row.addWidget(b)
            f.addRow("operator file", row)
            v.addLayout(f)
            root.addWidget(g)

            # advanced
            g = QGroupBox("Advanced")
            h = QHBoxLayout(g)
            box = QVBoxLayout()
            box.addWidget(QLabel("SigmaQ source"))
            self.c_sigmaq = QComboBox()
            self.c_sigmaq.addItems(["blend", "kfold", "kcal", "omit", "work"])
            box.addWidget(self.c_sigmaq)
            h.addLayout(box)
            box = QVBoxLayout()
            box.addWidget(QLabel("hist.lib"))
            self.e_histlib = QLineEdit()
            b = QPushButton("Browse…")
            b.clicked.connect(self.browse_histlib)
            row = QHBoxLayout()
            row.addWidget(self.e_histlib, 1)
            row.addWidget(b)
            box.addLayout(row)
            h.addLayout(box)
            self.cb_nohist = QCheckBox("no histogram library")
            h.addWidget(self.cb_nohist)
            h.addStretch(1)
            root.addWidget(g)

            # command preview
            self.l_cmd = QLabel("oasis-dm --mtzin … --mtzout …")
            self.l_cmd.setTextInteractionFlags(Qt.TextSelectableByMouse)
            self.l_cmd.setWordWrap(True)
            mono = QFont("Menlo")
            mono.setStyleHint(QFont.Monospace)
            self.l_cmd.setFont(mono)
            self.l_cmd.setStyleSheet(
                "background: #f4f6f9; border: 1px solid #d8dee6;"
                " border-radius: 6px; padding: 6px;")
            root.addWidget(self.l_cmd)

            # run row
            h = QHBoxLayout()
            self.b_run = QPushButton("Run")
            self.b_run.setStyleSheet(
                "background: #0369a1; color: white; font-weight: 600;"
                " padding: 6px 26px;")
            self.b_run.clicked.connect(self.run)
            self.b_stop = QPushButton("Stop")
            self.b_stop.setEnabled(False)
            self.b_stop.clicked.connect(self.stop)
            self.l_status = QLabel("ready")
            h.addWidget(self.b_run)
            h.addWidget(self.b_stop)
            h.addStretch(1)
            h.addWidget(self.l_status)
            root.addLayout(h)

            # chart + log
            self.chart = CcChart()
            self.chart.setMinimumHeight(140)
            root.addWidget(self.chart)
            self.log = QPlainTextEdit()
            self.log.setReadOnly(True)
            self.log.setFont(mono)
            self.log.setStyleSheet(
                "background: #10151c; color: #cfe3f2; font-family: Menlo;")
            root.addWidget(self.log, 1)

            # wire preview updates
            for w in (self.e_in, self.e_out, self.e_ncs, self.e_histlib):
                w.textChanged.connect(self.sync_cmd)
            for c in list(self.cols.values()) + [
                    self.c_mask, self.c_hist, self.c_metric, self.c_sigmaq]:
                c.currentTextChanged.connect(self.sync_cmd)
            for cb in (self.cb_stop, self.cb_extend, self.cb_aniso,
                       self.cb_pershell, self.cb_nohist):
                cb.toggled.connect(self.sync_cmd)
            for r in self.r_ncs.values():
                r.toggled.connect(self.sync_cmd)
            self.s_solc.valueChanged.connect(self.sync_cmd)
            self.s_cycles.valueChanged.connect(self.sync_cmd)
            self.cb_nohist.toggled.connect(self._sync_nohist)

        # -- input handling ------------------------------------------------

        def _debounce_input(self):
            # typing a path by hand: re-detect when the field settles
            path = self.e_in.text().strip()
            if hasattr(self, "_in_timer"):
                self._in_timer.stop()
                self._in_timer.deleteLater()
            from PySide6.QtCore import QTimer
            self._in_timer = QTimer(self)
            self._in_timer.setSingleShot(True)
            self._in_timer.timeout.connect(self._after_input)
            self._in_timer.start(400 if path else 0)

        def browse_in(self):
            path, _ = QFileDialog.getOpenFileName(
                self, "Choose input MTZ",
                self.e_in.text() or str(Path.home()),
                "MTZ files (*.mtz);;all files (*)")
            if not path:
                return
            self.e_in.setText(path)      # textChanged -> _after_input

        def _after_input(self):
            path = self.e_in.text().strip()
            if not path:
                self.l_info.setText("no input selected yet")
                return
            try:
                info = read_mtz_summary(path)
            except Exception as exc:                # noqa: BLE001
                self.l_info.setText("cannot read: %s" % exc)
                return
            for key, _label, wanted in COLUMNS:
                box = self.cols[key]
                box.clear()
                box.addItems(info["labels"] + ["(none)"])
                box.setCurrentText(pick_default(info["labels"], wanted, ""))
            self.l_info.setText(
                "%s · %s · %d reflections · d_min %.2f Å"
                % (Path(path).name, info["spacegroup"], info["nref"],
                   info["dmin"]))
            if not self.e_out.text().strip():
                p = Path(path)
                self.e_out.setText(str(p.with_name(p.stem + "_dm.mtz")))
            self.sync_cmd()

        def browse_out(self):
            path, _ = QFileDialog.getSaveFileName(
                self, "Choose output MTZ",
                self.e_out.text() or "dm.mtz", "MTZ files (*.mtz)")
            if path:
                self.e_out.setText(path)

        def browse_ncs(self):
            path, _ = QFileDialog.getOpenFileName(
                self, "Choose NCS operator file", "",
                "operator files (*.txt);;all files (*)")
            if path:
                self.r_ncs["file"].setChecked(True)
                self.e_ncs.setText(path)

        def browse_histlib(self):
            path, _ = QFileDialog.getOpenFileName(
                self, "Choose hist.lib", "", "all files (*)")
            if path:
                self.cb_nohist.setChecked(False)
                self.e_histlib.setText(path)

        # -- state ---------------------------------------------------------

        def state(self):
            ncs = next(k for k, r in self.r_ncs.items() if r.isChecked())
            st = {
                "mtzin": self.e_in.text().strip(),
                "mtzout": self.e_out.text().strip(),
                "solc": self.s_solc.value(),
                "cycles": self.s_cycles.value(),
                "auto_stop": self.cb_stop.isChecked(),
                "mask": self.c_mask.currentText(),
                "hist_mode": self.c_hist.currentText(),
                "free_metric": self.c_metric.currentText(),
                "extend": self.cb_extend.isChecked(),
                "aniso": self.cb_aniso.isChecked(),
                "per_shell_d": self.cb_pershell.isChecked(),
                "ncs": ncs,
                "ncs_file": self.e_ncs.text().strip() if ncs == "file" else "",
                "no_hist": self.cb_nohist.isChecked(),
                "hist_lib": "" if self.cb_nohist.isChecked()
                            else self.e_histlib.text().strip(),
                "sigmaq_source": self.c_sigmaq.currentText(),
            }
            for key, box in self.cols.items():
                st[key] = "" if box.currentText() == "(none)" \
                    else box.currentText()
            st["tag"] = (Path(st["mtzin"]).stem if st["mtzin"] else "dm")
            return st

        def _sync_ncs(self):
            mode = next((k for k, r in self.r_ncs.items()
                         if r.isChecked()), "none")
            self.e_ncs.setEnabled(mode == "file")

        def _sync_nohist(self, on):
            self.e_histlib.setDisabled(on)
            if on:
                self.e_histlib.clear()

        def sync_cmd(self):
            self.l_cmd.setText("$ " + command_string(self.state()))

        # -- run / stop ------------------------------------------------------

        def run(self):
            from PySide6.QtWidgets import QMessageBox
            st = self.state()
            problems = validate_state(st)
            if problems:
                QMessageBox.warning(self, "Cannot start",
                                    "· " + "\n· ".join(problems))
                return
            out = Path(st["mtzout"])
            if str(out.parent) and not out.parent.exists():
                out.parent.mkdir(parents=True, exist_ok=True)
            self.cycles = []
            self.chart.set_cycles([])
            self.log.clear()
            self.log.appendPlainText("$ " + command_string(st))
            self.b_run.setEnabled(False)
            self.b_stop.setEnabled(True)
            self.l_status.setText("running…")
            self.l_status.setStyleSheet("")

            from PySide6.QtCore import QProcess, QProcessEnvironment
            self.proc = QProcess(self)
            self.proc.setProcessChannelMode(QProcess.MergedChannels)
            env = QProcessEnvironment.systemEnvironment()
            ppath = _pythonpath_env().get("PYTHONPATH", "")
            if ppath:
                env.insert("PYTHONPATH", ppath)
            self.proc.setProcessEnvironment(env)
            self.proc.readyReadStandardOutput.connect(self._on_output)
            self.proc.finished.connect(self._on_finished)
            self._t0 = time.time()
            self.proc.start(sys.executable,
                            ["-m", "oasis_dm"] + build_argv(st))

        def _on_output(self):
            while True:
                line = self.proc.readLine()
                if not line:
                    break
                text = bytes(line).decode("utf-8", "replace").rstrip("\n")
                self.log.appendPlainText(text)
                hit = parse_cycle_line(text)
                if hit and (not self.cycles or hit[0] > self.cycles[-1][0]):
                    self.cycles.append(hit)
                    self.chart.set_cycles(self.cycles)
                if self.cycles:
                    self.l_status.setText(
                        "running… cycle %d · %d s"
                        % (self.cycles[-1][0], int(time.time() - self._t0)))

        def _on_finished(self, rc, _status):
            self.b_run.setEnabled(True)
            self.b_stop.setEnabled(False)
            if rc == 0:
                self.l_status.setText("done in %d s"
                                       % int(time.time() - self._t0))
                self.l_status.setStyleSheet(
                    "color: #059669; font-weight: 600;")
            else:
                self.l_status.setText("FAILED (exit %d) — see the log" % rc)
                self.l_status.setStyleSheet(
                    "color: #b91c1c; font-weight: 600;")

        def stop(self):
            from PySide6.QtCore import QProcess
            if self.proc is not None \
                    and self.proc.state() != QProcess.NotRunning:
                self.proc.kill()
                self.log.appendPlainText("[gui] stopped by user")

        def closeEvent(self, event):
            self.stop()
            super().closeEvent(event)

    return Win()


def main(argv=None):                                 # noqa: ARG001 (CLI parity)
    if not qt_available():
        sys.stderr.write(
            "oasis-dm-gui needs PySide6 (the desktop Qt bindings).\n"
            "Install it with:  pip install PySide6\n"
            "or reinstall the package with:  pip install \".[gui]\"\n")
        return 1
    app, win = build_app()
    win.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
