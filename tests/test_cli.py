"""Smoke tests: the CLI runs, the column contract holds, and hist.lib is optional.

These tests deliberately avoid a hard dependency on pytest so they can run in a
minimal offline environment:

    python tools/run_tests.py          # no pytest required
    pytest tests/                      # if pytest is available

They use only plain asserts, `tempfile` and a local `raises` helper, so pytest
collects them unchanged.
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ENV = {**os.environ, "PYTHONPATH": str(ROOT / "src")}


class raises:
    """Minimal stand-in for pytest.raises (context manager)."""

    def __init__(self, expected, match=None):
        self.expected = expected
        self.match = match

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        if exc_type is None:
            raise AssertionError("expected %s, nothing was raised" % self.expected)
        if not issubclass(exc_type, self.expected):
            return False
        if self.match is not None and self.match not in str(exc):
            raise AssertionError("message %r does not contain %r" % (str(exc), self.match))
        return True


def run(*args):
    return subprocess.run([sys.executable, "-m", "oasis_dm", *args],
                          cwd=str(ROOT), capture_output=True, text=True, env=ENV)


def test_version_flag():
    result = run("--version")
    assert result.returncode == 0, result.stderr
    assert "oasis-dm" in result.stdout


def test_runs_without_a_histogram_library():
    import gemmi

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        out = tmp / "out.mtz"
        result = run("--mtzin", "tests/data/3o2e.mtz", "--mtzout", str(out),
                     "--stats", str(tmp / "stats.json"), "--no-hist",
                     "--no-auto-stop", "--cycles", "2", "--quiet")
        assert result.returncode == 0, result.stderr
        labels = set(gemmi.read_mtz_file(str(out)).column_labels())
        assert {"DM_FDM", "DM_PHIDM", "DM_FOMDM", "DM_FWT", "DM_PHWT"} <= labels
        assert (tmp / "stats.json").is_file()


def test_missing_hist_lib_is_reported_clearly():
    sys.path.insert(0, str(ROOT / "src"))
    import oasis_dm.paths as paths

    saved = {k: os.environ.pop(k, None) for k in ("OASIS_DM_HIST_LIB", "CLIBD")}
    cwd = os.getcwd()
    real_bundled = paths.bundled_path
    try:
        with tempfile.TemporaryDirectory() as tmp:
            os.chdir(tmp)
            # the package ships a fallback library; pretend it is absent to test
            # the "nothing found" path
            paths.bundled_path = lambda: Path(tmp) / "nope.lib"
            with raises(ValueError, match="hist.lib"):
                paths.resolve_hist_lib()
            with raises(ValueError, match="not a file"):
                paths.resolve_hist_lib("/nonexistent/hist.lib")
            # and with the fallback in place it resolves without any CCP4 install
            paths.bundled_path = real_bundled
            assert Path(paths.resolve_hist_lib()).is_file()
    finally:
        paths.bundled_path = real_bundled
        os.chdir(cwd)
        for key, value in saved.items():
            if value is not None:
                os.environ[key] = value


def test_legacy_handdm_script_name_still_works():
    """OASIS deployments resolve the engine as a script named handdm_v2.py."""
    result = subprocess.run([sys.executable, str(ROOT / "handdm_v2.py"), "--version"],
                            cwd=str(ROOT), capture_output=True, text=True, env=ENV)
    assert result.returncode == 0, result.stderr
    assert "oasis-dm" in result.stdout


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
            print("FAIL %s: %s" % (name, exc))
    print("\n%s" % ("all tests passed" if not failures else "%d test(s) failed" % failures))
    sys.exit(1 if failures else 0)
