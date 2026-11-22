"""Browser-based GUI for oasis_dm.

``oasis-dm-gui`` (or ``python -m oasis_dm.gui``) starts a small local HTTP
server (standard library only, bound to 127.0.0.1) and opens the default
browser.  The engine itself always runs as a *subprocess* of
``python -m oasis_dm``, so a run can be stopped cleanly and its log
streams into the page.

Why not Tkinter?  It would add no dependency on paper, but in practice the
Tcl/Tk shipped with system Pythons is frequently broken (on macOS 15 with
the Xcode toolchain Python it aborts the process outright), and the page
needs no GUI toolkit on the Python side at all.  Everything the browser
needs is served from `gui_page.py`; there are no external assets and no
JavaScript dependencies.

All POST requests must carry the per-session token (in ``X-Oasis-Token``
and in the URL path), which keeps other local processes and cross-origin
pages from driving the engine.
"""
from __future__ import annotations

import json
import os
import re
import shlex
import subprocess
import sys
import threading
import time
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

from . import __version__
from .gui_page import PAGE

# --- pure logic (unit-tested in tests/test_gui.py) -------------------------


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
    """Subprocess env that finds oasis_dm from a source checkout."""
    env = dict(os.environ)
    here = Path(__file__).resolve().parent          # .../src/oasis_dm
    src = here.parent                                # .../src
    if (src / "oasis_dm" / "__init__.py").exists() and str(src) \
            not in env.get("PYTHONPATH", ""):
        env["PYTHONPATH"] = str(src) + os.pathsep + env.get("PYTHONPATH", "")
    return env


def list_dir(path):
    """Entries for the in-page file browser (directories first)."""
    p = Path(path).expanduser().resolve()
    if not p.is_dir():
        raise FileNotFoundError(path)
    entries = []
    for it in sorted(p.iterdir(), key=lambda e: e.name.lower()):
        if it.name.startswith("."):
            continue
        is_dir = it.is_dir()
        entries.append({
            "name": it.name,
            "dir": is_dir,
            "mtz": (not is_dir) and it.suffix.lower() == ".mtz",
        })
    return {
        "cwd": str(p),
        "parent": str(p.parent) if p.parent != p else None,
        "entries": entries,
    }


# --- job management ---------------------------------------------------------


class Job:
    """One engine subprocess plus its captured log."""

    def __init__(self, argv):
        self.argv = argv
        self.proc = subprocess.Popen(
            [sys.executable, "-m", "oasis_dm"] + argv,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
            env=_pythonpath_env())
        self.lines = []
        self.cycles = []                 # [[cyc, meanFOM, freeCC], ...]
        self.rc = None
        self.t0 = time.time()
        threading.Thread(target=self._reader, daemon=True).start()

    def _reader(self):
        try:
            for line in self.proc.stdout:
                self.lines.append(line)
                hit = parse_cycle_line(line)
                if hit and (not self.cycles or hit[0] > self.cycles[-1][0]):
                    self.cycles.append(hit)
        finally:
            self.rc = self.proc.wait()

    def poll(self, seq):
        return {
            "running": self.rc is None,
            "rc": self.rc,
            "total": len(self.lines),
            "lines": self.lines[seq:],
            "cycles": self.cycles,
            "elapsed": round(time.time() - self.t0, 1),
        }

    def stop(self):
        if self.proc.poll() is None:
            self.proc.terminate()
            self.lines.append("\n[gui] stopped by user\n")


class GuiState:
    def __init__(self):
        self.job = None
        self.mtzout = None
        self.lock = threading.Lock()

    def start(self, argv, mtzout):
        with self.lock:
            if self.job is not None and self.job.rc is None:
                return False
            self.job = Job(argv)
            self.mtzout = mtzout
            return True

    def status(self):
        with self.lock:
            job = self.job
            if job is None:
                return None
            info = job.poll(10 ** 9)      # everything but the lines
            info["lines"] = []
            info["output_exists"] = bool(
                self.mtzout and Path(self.mtzout).is_file())
            return info


def make_handler(state, token):
    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *a):        # noqa: ARG002 (quiet by default)
            pass

        # -- plumbing ----------------------------------------------------

        def _authorized(self):
            ok = (self.headers.get("X-Oasis-Token") == token
                  or self.path.startswith("/t/%s/" % token))
            if not ok:
                self._json(403, {"ok": False, "error": "bad token"})
            return ok

        def _json(self, code, payload):
            body = json.dumps(payload).encode()
            self.send_response(code)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(body)

        def _body(self):
            n = int(self.headers.get("Content-Length") or 0)
            return json.loads(self.rfile.read(n) or b"{}")

        # -- routes ------------------------------------------------------

        def do_GET(self):                 # noqa: N802 (http.server API)
            if self.path in ("/", "/index.html") \
                    or re.fullmatch(
                        r"/t/%s(/|/index\.html)?" % re.escape(token),
                        self.path):
                body = PAGE.replace("__TOKEN__", token).encode()
                self.send_response(200)
                self.send_header("Content-Type", "text/html; charset=utf-8")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
                return
            if not self._authorized():
                return
            if self.path.startswith("/api/poll"):
                info = state.status()
                if info is None:
                    self._json(200, {"running": False, "idle": True})
                    return
                seq = 0
                q = self.path.partition("?")[2]
                for kv in q.split("&"):
                    if kv.startswith("seq="):
                        seq = int(kv.partition("=")[2] or 0)
                info = state.job.poll(seq)
                info["idle"] = False
                info["output_exists"] = bool(
                    state.mtzout and Path(state.mtzout).is_file())
                self._json(200, info)
                return
            if self.path.startswith("/api/list"):
                from urllib.parse import parse_qs
                q = parse_qs(self.path.partition("?")[2])
                path = (q.get("path") or [os.getcwd()])[0]
                try:
                    self._json(200, {"ok": True, **list_dir(path)})
                except (FileNotFoundError, PermissionError) as exc:
                    self._json(200, {"ok": False, "error": str(exc)})
                return
            self._json(404, {"ok": False, "error": "not found"})

        def do_POST(self):                # noqa: N802 (http.server API)
            if not self._authorized():
                return
            try:
                body = self._body()
            except (ValueError, json.JSONDecodeError):
                self._json(400, {"ok": False, "error": "bad JSON"})
                return
            if self.path.startswith("/api/open"):
                try:
                    info = read_mtz_summary(body["path"])
                except Exception as exc:                # noqa: BLE001
                    self._json(200, {"ok": False, "error": str(exc)})
                    return
                self._json(200, {"ok": True, **info})
                return
            if self.path.startswith("/api/list"):
                try:
                    self._json(200, {"ok": True,
                                     **list_dir(body.get("path") or ".")})
                except (FileNotFoundError, PermissionError) as exc:
                    self._json(200, {"ok": False, "error": str(exc)})
                return
            if self.path.startswith("/api/preview"):
                st = body.get("state", {})
                self._json(200, {"ok": True,
                                 "cmd": command_string(st)})
                return
            if self.path.startswith("/api/run"):
                st = body.get("state", {})
                problems = validate_state(st)
                if problems:
                    self._json(200, {"ok": False, "problems": problems})
                    return
                out = Path(st["mtzout"])
                if str(out.parent) and not out.parent.exists():
                    out.parent.mkdir(parents=True, exist_ok=True)
                argv = build_argv(st)
                if not state.start(argv, str(out)):
                    self._json(200, {"ok": False,
                                     "problems": ["a run is already active"]})
                    return
                self._json(200, {"ok": True, "cmd": command_string(st)})
                return
            if self.path.startswith("/api/stop"):
                with state.lock:
                    if state.job is not None:
                        state.job.stop()
                self._json(200, {"ok": True})
                return
            self._json(404, {"ok": False, "error": "not found"})

    return Handler


def command_string(st):
    return "oasis-dm " + " ".join(shlex.quote(a) for a in build_argv(st))


def serve(open_browser=True, host="127.0.0.1", port=0):
    """Start the GUI server (blocks).  Returns only after the server exits."""
    import secrets

    token = secrets.token_urlsafe(18)
    state = GuiState()
    httpd = ThreadingHTTPServer((host, port),
                                make_handler(state, token))
    host_, port_ = httpd.server_address[:2]
    url = "http://%s:%d/t/%s/" % (host_, port_, token)
    print("oasis-dm %s GUI: %s  (Ctrl-C to quit)" % (__version__, url))
    if open_browser:
        threading.Timer(0.3, lambda: webbrowser.open(url)).start()
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        with state.lock:
            if state.job is not None:
                state.job.stop()
        httpd.server_close()
    return 0


def main(argv=None):                                 # noqa: ARG001
    return serve()


if __name__ == "__main__":
    sys.exit(main())
