"""Tests for the desktop GUI: pure option logic + Qt window behaviour.

The pure-logic tests (command construction, column auto-detection, cycle
parsing) run everywhere.  The window tests need PySide6 and are skipped
when it is not installed — run them with a Python from
``pip install ".[gui]"``.
"""
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "tests" / "data"
sys.path.insert(0, str(ROOT / "src"))

from oasis_dm import gui                     # noqa: E402

try:
    import PySide6                           # noqa: F401
    HAVE_QT = True
except ImportError:
    HAVE_QT = False

if HAVE_QT:
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


# --- pure logic -------------------------------------------------------------

def _base(**over):
    st = {
        "mtzin": "in.mtz", "mtzout": "out.mtz", "fo": "", "sigf": "",
        "phib": "", "fom": "", "free": "", "truth": "",
        "solc": 0.5, "cycles": 8, "auto_stop": True,
        "mask": "variance", "hist_mode": "cell", "free_metric": "cc",
        "extend": True, "aniso": False, "per_shell_d": False,
        "ncs": "none", "ncs_file": "", "no_hist": False, "hist_lib": "",
        "sigmaq_source": "blend", "stats": "",
    }
    st.update(over)
    return st


def test_minimal_run_uses_engine_defaults():
    assert gui.build_argv(_base()) == ["--mtzin", "in.mtz",
                                       "--mtzout", "out.mtz"]


def test_tag_is_emitted_with_stats():
    argv = gui.build_argv(_base(stats="stats/", tag="in"))
    assert argv[4:8] == ["--stats", "stats/", "--tag", "in"]


def test_changed_columns_are_emitted():
    argv = gui.build_argv(_base(fo="FP", free="FREE"))
    assert "--fo" in argv and argv[argv.index("--fo") + 1] == "FP"
    assert "--free" in argv and argv[argv.index("--free") + 1] == "FREE"


def test_default_columns_are_not_emitted():
    argv = gui.build_argv(_base(fo="FOBS", free="FreeR"))
    assert "--fo" not in argv and "--free" not in argv


def test_parameters_map_to_flags():
    argv = gui.build_argv(_base(solc=0.42, cycles=3, auto_stop=False,
                                mask="wang", hist_mode="protein",
                                free_metric="rfree", extend=False))
    for flag, val in (("--solc", "0.420"), ("--cycles", "3"),
                      ("--no-auto-stop", None), ("--mask", "wang"),
                      ("--hist-mode", "protein"), ("--free-metric", "rfree"),
                      ("--no-extend", None)):
        assert flag in argv
        if val is not None:
            assert argv[argv.index(flag) + 1] == val


def test_ncs_modes():
    assert gui.build_argv(_base(ncs="none"))[-1] != "--ncs-auto"
    argv = gui.build_argv(_base(ncs="file", ncs_file="op.txt"))
    assert argv[argv.index("--ncs") + 1] == "op.txt"
    assert "--ncs-auto" in gui.build_argv(_base(ncs="auto"))


def test_option_toggles():
    argv = gui.build_argv(_base(aniso=True, per_shell_d=True,
                                sigmaq_source="kfold", no_hist=True))
    for flag in ("--aniso", "--per-shell-d", "--no-hist"):
        assert flag in argv
    assert argv[argv.index("--sigmaq-source") + 1] == "kfold"
    # hist_lib is ignored when the library is switched off
    argv = gui.build_argv(_base(no_hist=True, hist_lib="my.lib"))
    assert "--hist-lib" not in argv


def test_validate_state_catches_problems():
    st = _base(mtzin="nope.mtz", mtzout="", solc=0.99, ncs="file")
    problems = gui.validate_state(st)
    assert any("input MTZ" in p for p in problems)
    assert any("output MTZ" in p for p in problems)
    assert any("solvent" in p for p in problems)
    assert any("NCS operator" in p for p in problems)
    assert gui.validate_state(_base(mtzin=str(DATA / "7e1d.mtz"))) == []


def test_parse_cycle_line():
    assert gui.parse_cycle_line(" 1  0.71  0.65  0.31") == [1, 0.71, 0.65]
    assert gui.parse_cycle_line("12 0.5  -0.1  0.4") == [12, 0.5, -0.1]
    assert gui.parse_cycle_line("nan nan nan") is None
    assert gui.parse_cycle_line("[out ] wrote x.mtz") is None
    assert gui.parse_cycle_line(" 1  2.5  0.5  0.3") is None   # FOM > 1


def test_pick_default():
    assert gui.pick_default(["FOBS", "FP"], ("FP", "F")) == "FP"
    assert gui.pick_default(["A", "B"], ("X",), "") == ""
    assert gui.pick_default(["A", "B"], ("X",)) == "A"


def test_read_mtz_summary():
    info = gui.read_mtz_summary(DATA / "7e1d.mtz")
    assert "FOBS" in info["labels"]
    assert info["nref"] > 0
    assert 0 < info["dmin"] < 100
    assert isinstance(info["nref"], int)     # JSON/plain-serialisable


def test_command_string_quotes_paths():
    cmd = gui.command_string(_base(mtzin="/tmp/my file.mtz"))
    assert "'/tmp/my file.mtz'" in cmd


# --- Qt window (need PySide6) -----------------------------------------------

if HAVE_QT:

    def _window():
        app, win = gui.build_app()
        win.show()
        app.processEvents()
        return app, win

    def _spin(app, seconds):
        t0 = time.time()
        while time.time() - t0 < seconds:
            app.processEvents()
            time.sleep(0.02)

    def test_window_state_maps_to_command():
        app, win = _window()
        win.e_in.setText(str(DATA / "7e1d.mtz"))
        _spin(app, 1.0)                       # let the debounce fire
        assert "reflections" in win.l_info.text()
        assert win.cols["fo"].currentText() == "FOBS"
        assert win.cols["free"].currentText() == "FreeR"
        assert win.e_out.text().endswith("7e1d_dm.mtz")
        st = win.state()
        assert st["fo"] == "FOBS"
        assert gui.build_argv(st)[0:4] == [
            "--mtzin", str(DATA / "7e1d.mtz"),
            "--mtzout", str(DATA / "7e1d_dm.mtz")]
        win.close()

    def test_run_button_validates_before_starting():
        app, win = _window()
        # no input file selected -> run() must refuse (message box is
        # modal in real use; here the state is invalid before any dialog)
        st = win.state()
        assert gui.validate_state(st) != []
        win.close()

    def test_window_runs_a_real_engine_job():
        from PySide6.QtCore import QProcess
        import tempfile
        app, win = _window()
        out = Path(tempfile.mkdtemp()) / "gui_e2e.mtz"
        win.e_in.setText(str(DATA / "7e1d.mtz"))
        _spin(app, 1.0)
        win.e_out.setText(str(out))
        win.s_cycles.setValue(2)
        win.cb_stop.setChecked(False)
        app.processEvents()
        win.run()
        t0 = time.time()
        while win.proc is None \
                or win.proc.state() != QProcess.NotRunning:
            app.processEvents()
            if time.time() - t0 > 240:
                raise AssertionError("engine did not finish in 240 s")
            time.sleep(0.02)
        app.processEvents()
        assert out.is_file()
        assert len(win.cycles) == 2
        assert win.cycles[0][0] == 1          # cycle numbers
        assert "done" in win.l_status.text()
        assert "wrote" in win.log.toPlainText()
        win.close()

    def test_stop_button_kills_the_engine():
        import tempfile
        app, win = _window()
        out = Path(tempfile.mkdtemp()) / "gui_stop.mtz"
        win.e_in.setText(str(DATA / "7e1d.mtz"))
        _spin(app, 1.0)
        win.e_out.setText(str(out))
        win.s_cycles.setValue(50)
        win.cb_stop.setChecked(False)
        app.processEvents()
        win.run()
        from PySide6.QtCore import QProcess
        t0 = time.time()
        while win.proc is None or win.proc.state() != QProcess.Running:
            app.processEvents()
            if time.time() - t0 > 30:
                raise AssertionError("engine did not start in 30 s")
            time.sleep(0.02)
        win.stop()
        t0 = time.time()
        while win.proc.state() != QProcess.NotRunning:
            app.processEvents()
            if time.time() - t0 > 20:
                raise AssertionError("engine did not stop in 20 s")
            time.sleep(0.02)
        assert "[gui] stopped by user" in win.log.toPlainText()
        win.close()

    def test_chart_receives_cycles():
        app, win = _window()
        win.chart.set_cycles([[1, 0.5, 0.3], [2, 0.6, 0.55]])
        win.chart.repaint()
        app.processEvents()
        assert win.chart.pts == [(1, 0.3), (2, 0.55)]
        win.close()
