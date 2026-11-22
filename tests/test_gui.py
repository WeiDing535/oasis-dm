"""GUI logic and HTTP end-to-end tests.

The argv construction, column auto-detection, log parsing, directory
listing and the local HTTP server are all exercised here.  The final test
starts a real server in a thread and drives a real engine run through the
same JSON API the browser page uses.
"""
import json
import sys
import tempfile
import threading
import time
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from oasis_dm.gui import (          # noqa: E402
    GuiState, build_argv, command_string, list_dir, make_handler,
    parse_cycle_line, pick_default, read_mtz_summary, validate_state)


def _base(**over):
    st = dict(
        mtzin="/tmp/in.mtz", mtzout="/tmp/out.mtz",
        fo="FOBS", sigf="SIGFOBS", phib="PHIB", fom="FOM",
        free="FreeR", truth="PHIC",
        solc=0.5, cycles=8, auto_stop=True,
        mask="variance", hist_mode="cell", free_metric="cc",
        extend=True, aniso=False, per_shell_d=False,
        ncs="none", ncs_file="", no_hist=False, hist_lib="",
        sigmaq_source="blend",
    )
    st.update(over)
    return st


def test_default_state_is_a_minimal_command():
    """Nothing changed from the defaults -> just files, no option spam."""
    argv = build_argv(_base())
    assert argv[:4] == ["--mtzin", "/tmp/in.mtz", "--mtzout", "/tmp/out.mtz"]
    assert "--solc" not in argv and "--ncs" not in argv
    assert "--tag" not in argv          # only meaningful with --stats


def test_tag_is_emitted_with_stats():
    argv = build_argv(_base(stats="/tmp/s.json", tag="myset"))
    assert argv[argv.index("--stats") + 1] == "/tmp/s.json"
    assert argv[argv.index("--tag") + 1] == "myset"


def test_changed_columns_and_parameters_are_emitted():
    argv = build_argv(_base(fo="FP", free="FREE", solc=0.42, cycles=12,
                            auto_stop=False, mask="wang"))
    for flag, val in (("--fo", "FP"), ("--free", "FREE"),
                      ("--solc", "0.420"), ("--cycles", "12"),
                      ("--mask", "wang")):
        assert argv[argv.index(flag) + 1] == val, (flag, argv)
    assert "--no-auto-stop" in argv
    # unchanged column defaults are not spelled out
    assert "--phib" not in argv


def test_none_column_is_omitted():
    argv = build_argv(_base(truth=""))
    assert "--truth" not in argv


def test_ncs_modes_map_to_flags():
    assert "--ncs" not in build_argv(_base())
    argv = build_argv(_base(ncs="file", ncs_file="/tmp/op.txt"))
    assert argv[argv.index("--ncs") + 1] == "/tmp/op.txt"
    assert "--ncs-auto" not in argv
    argv = build_argv(_base(ncs="auto"))
    assert "--ncs-auto" in argv and "--ncs" not in argv


def test_option_toggles():
    argv = build_argv(_base(extend=False, aniso=True, per_shell_d=True,
                            no_hist=True, sigmaq_source="kfold"))
    for flag in ("--no-extend", "--aniso", "--per-shell-d", "--no-hist"):
        assert flag in argv, (flag, argv)
    assert argv[argv.index("--sigmaq-source") + 1] == "kfold"


def test_validate_state_catches_problems():
    assert validate_state(_base(mtzin="/no/such/file.mtz"))
    assert validate_state(_base(mtzout=""))
    assert validate_state(_base(solc=0.99))
    assert validate_state(_base(ncs="file", ncs_file=""))
    good = _base(mtzin=str(ROOT / "tests" / "data" / "7e1d.mtz"))
    assert validate_state(good) == []


def test_parse_cycle_line():
    """The engine's statistics table feeds the free-cc chart."""
    line = "  1  0.7123  0.6544  0.2810   15203.94   0.147   20.49  0.8912"
    assert parse_cycle_line(line) == [1, 0.7123, 0.6544]
    line8 = ("  8  0.6698  0.7182  0.6046  15231.38   0.143  17.79   0.9152"
             "      *")
    assert parse_cycle_line(line8) == [8, 0.6698, 0.7182]
    # header and message lines must not parse
    assert parse_cycle_line("cyc  meanFOM  freeCC  freeR ...") is None
    assert parse_cycle_line("[v2  ] NCS gate: free cc 0.786 ...") is None
    # NaN columns (no free set) are rejected
    assert parse_cycle_line("  1   nan    nan    nan ...") is None


def test_pick_default_prefers_common_labels():
    assert pick_default(["FP", "SIGFP", "PHIB"], ["FOBS", "FP"]) == "FP"
    assert pick_default(["A", "B"], ["FOBS"]) == "A"
    assert pick_default(["A", "B"], ["FOBS"], fallback="") == ""


def test_read_mtz_summary_reports_header():
    info = read_mtz_summary(ROOT / "tests" / "data" / "7e1d.mtz")
    assert info["nref"] > 0
    assert info["dmin"] < 3.0
    assert "FOBS" in info["labels"] and "FreeR" in info["labels"]


def test_list_dir_sorts_and_marks_mtz():
    with tempfile.TemporaryDirectory() as td:
        (Path(td) / "b.mtz").write_bytes(b"")
        (Path(td) / "a.txt").write_bytes(b"")
        (Path(td) / "sub").mkdir()
        info = list_dir(td)
        names = [e["name"] for e in info["entries"]]
        assert names == ["a.txt", "b.mtz", "sub"]   # case-insensitive sort
        assert info["entries"][1]["mtz"] and not info["entries"][0]["mtz"]
        assert info["parent"]          # tmp has a parent


def _serve(state, token):
    from http.server import ThreadingHTTPServer
    httpd = ThreadingHTTPServer(("127.0.0.1", 0), make_handler(state, token))
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    return httpd


def _req(port, method, path, token, body=None):
    data = json.dumps(body).encode() if body is not None else None
    r = urllib.request.Request(
        "http://127.0.0.1:%d%s" % (port, path), data=data, method=method,
        headers={"X-Oasis-Token": token,
                 "Content-Type": "application/json"})
    with urllib.request.urlopen(r, timeout=30) as resp:
        return json.loads(resp.read())


def test_http_api_rejects_bad_token():
    state = GuiState()
    httpd = _serve(state, "sekrit")
    try:
        try:
            _req(httpd.server_address[1], "GET", "/api/poll", "wrong")
            ok = False
        except urllib.error.HTTPError as exc:
            ok = exc.code == 403
        assert ok
    finally:
        httpd.shutdown()


def test_http_api_open_preview_and_list():
    state = GuiState()
    httpd = _serve(state, "sekrit")
    port = httpd.server_address[1]
    try:
        r = _req(port, "POST", "/api/open", "sekrit",
                 {"path": str(ROOT / "tests" / "data" / "7e1d.mtz")})
        assert r["ok"] and "FOBS" in r["labels"] and r["nref"] > 0
        r = _req(port, "POST", "/api/open", "sekrit",
                 {"path": str(ROOT / "nope.mtz")})
        assert not r["ok"]
        r = _req(port, "POST", "/api/preview", "sekrit",
                 {"state": _base(cycles=3, mask="wang")})
        assert r["ok"] and "--cycles 3" in r["cmd"] and "--mask wang" \
            in r["cmd"]
        assert "$" not in r["cmd"]       # no shell decoration server-side
        r = _req(port, "POST", "/api/list", "sekrit",
                 {"path": str(ROOT / "tests" / "data")})
        assert r["ok"] and any(e["mtz"] for e in r["entries"])
        # and the page itself is served with the token embedded
        with urllib.request.urlopen(
                "http://127.0.0.1:%d/" % port, timeout=30) as resp:
            page = resp.read().decode()
        assert "oasis-dm" in page and "sekrit" in page \
            and "__TOKEN__" not in page
    finally:
        httpd.shutdown()


def test_run_validation_via_http():
    state = GuiState()
    httpd = _serve(state, "sekrit")
    try:
        r = _req(httpd.server_address[1], "POST", "/api/run", "sekrit",
                 {"state": _base(mtzin="/no/such.mtz")})
        assert not r["ok"] and r["problems"]
    finally:
        httpd.shutdown()


def test_http_api_runs_a_real_engine_job():
    """Full round trip: run, poll the streaming log, check the output."""
    with tempfile.TemporaryDirectory() as td:
        out = str(Path(td) / "gui_e2e.mtz")
        st = _base(mtzin=str(ROOT / "tests" / "data" / "7e1d.mtz"),
                   mtzout=out, cycles=2, auto_stop=False)
        state = GuiState()
        httpd = _serve(state, "sekrit")
        port = httpd.server_address[1]
        try:
            r = _req(port, "POST", "/api/run", "sekrit", {"state": st})
            assert r["ok"], r
            # refuse a second concurrent run
            r2 = _req(port, "POST", "/api/run", "sekrit", {"state": st})
            assert not r2["ok"] and "already active" in r2["problems"][0]

            seq, saw_cycle, deadline = 0, False, time.time() + 600
            while time.time() < deadline:
                info = _req(port, "GET", "/api/poll?seq=%d" % seq, "sekrit")
                seq = info["total"]
                if info["cycles"]:
                    saw_cycle = True
                    cyc0 = info["cycles"][0]
                    assert cyc0[0] >= 1 and 0 < cyc0[1] <= 1
                if not info["running"]:
                    break
                time.sleep(1.0)
            assert not info["running"] and info["rc"] == 0, info
            assert saw_cycle and len(info["cycles"]) == 2
            assert info["output_exists"] and Path(out).is_file()
            assert any("wrote" in ln for ln in state.job.lines)
        finally:
            httpd.shutdown()
            if state.job is not None and state.job.rc is None:
                state.job.stop()


def test_http_api_stop_kills_a_job():
    """Stop must terminate a long run within a few seconds."""
    with tempfile.TemporaryDirectory() as td:
        st = _base(mtzin=str(ROOT / "tests" / "data" / "7e1d.mtz"),
                   mtzout=str(Path(td) / "s.mtz"), cycles=50,
                   auto_stop=False)
        state = GuiState()
        httpd = _serve(state, "sekrit")
        port = httpd.server_address[1]
        try:
            assert _req(port, "POST", "/api/run", "sekrit",
                        {"state": st})["ok"]
            time.sleep(3.0)               # let it get going
            assert state.job.rc is None   # still running
            assert _req(port, "POST", "/api/stop", "sekrit")["ok"]
            t0 = time.time()
            while state.job.rc is None and time.time() - t0 < 20:
                time.sleep(0.2)
            assert state.job.rc is not None
        finally:
            httpd.shutdown()
            if state.job is not None and state.job.rc is None:
                state.job.stop()


def test_command_string_is_quoted():
    s = command_string(_base(mtzin="/tmp/my file.mtz"))
    assert s.startswith("oasis-dm --mtzin ") and "'/tmp/my file.mtz'" in s


if __name__ == "__main__":
    failures = 0
    for name, fn in sorted(globals().items()):
        if not name.startswith("test_") or not callable(fn):
            continue
        try:
            fn()
            print("PASS %s" % name)
        except Exception as exc:  # noqa: BLE001
            failures += 1
            print("FAIL %s: %r" % (name, exc))
    sys.exit(1 if failures else 0)
