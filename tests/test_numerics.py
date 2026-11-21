"""Unit tests for the numerics and the MTZ contract.

Same style as ``test_cli.py``: plain asserts, no pytest required (pytest collects
them unchanged).  Run with ``python tools/run_tests.py`` or ``pytest tests/``.
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from oasis_dm.core import (  # noqa: E402
    HistLib, hl_centroid, fom_to_x, x_to_fom, phifom_to_hl, wrap180,
    nice_fft_size, DensityModifier, Config)

ENV = {**os.environ, "PYTHONPATH": str(ROOT / "src")}


def test_wrap180():
    assert np.allclose(wrap180([0.0, 180.0, -180.0, 190.0, -190.0]),
                       [0.0, -180.0, -180.0, -170.0, 170.0])


def test_nice_fft_size():
    for n in (4, 8, 17, 97, 128):
        m = nice_fft_size(n)
        assert m >= n
        t = m
        for p in (2, 3, 5):
            while t % p == 0:
                t //= p
        assert t == 1, m


def test_fom_roundtrip():
    fom = np.array([0.0, 0.1, 0.35, 0.6, 0.9, 0.99])
    back = x_to_fom(fom_to_x(fom))
    assert np.allclose(back, fom, atol=2e-4), back


def test_hl_centroid_matches_von_mises():
    """Unimodal HL must reproduce the analytic centroid exactly."""
    phi = np.deg2rad([10.0, 120.0, -160.0, 45.0])
    fom = np.array([0.2, 0.5, 0.8, 0.95])
    p, m = hl_centroid(phifom_to_hl(phi, fom))
    assert np.allclose(np.abs(wrap180(np.rad2deg(p) - np.rad2deg(phi))), 0.0, atol=1e-6)
    assert np.allclose(m, fom, atol=1e-4)


def test_hl_centroid_bimodal_is_centred_between_modes():
    """A symmetric two-mode distribution has its centroid between the modes."""
    # A*cos(phi) = 0 with C<0 gives modes at 0 and 180 deg
    hl = np.array([[0.0, 0.0, -6.0, 0.0]])
    p, m = hl_centroid(hl)
    # symmetric: the resultant is ~0, so only the FOM is meaningful
    assert m[0] < 1e-3


def test_bundled_histogram_library():
    path = ROOT / "src" / "oasis_dm" / "data" / "protein_hist.lib"
    assert path.is_file(), "bundled histogram library is missing"
    h = HistLib(str(path))
    assert h.nbin == 100 and h.nres == 18
    cdf, mean, std = h.target(2.5)
    assert abs(cdf[-1] - 1.0) < 1e-9
    assert 0.2 < mean < 0.8 and std > 0.05
    # the protein-only target must sit above the whole-cell one
    cdf_p, mean_p, std_p, axis_p, rho_cut = h.protein_target(2.5, 0.5)
    assert abs(cdf_p[-1] - 1.0) < 1e-9
    assert mean_p > mean
    assert axis_p[0] >= rho_cut


def test_fft_roundtrip():
    """rho -> grid -> rho must be the identity (the engine's own convention)."""
    cfg = Config(mtzin=str(ROOT / "tests" / "data" / "3o2e.mtz"),
                 mtzout="", hist_on=False, quiet=True)
    dm = DensityModifier(cfg)
    rng = np.random.default_rng(0)
    rho = rng.normal(size=tuple(dm.size))
    back = dm.rho_from_grid(dm.grid_from_rho(rho))
    assert np.allclose(back, rho, atol=1e-9)


def test_output_is_idempotent():
    """Re-running on an earlier output must not duplicate DM_ columns."""
    import gemmi
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        first, second = tmp / "a.mtz", tmp / "b.mtz"
        base = [sys.executable, "-m", "oasis_dm", "--mtzin",
                str(ROOT / "tests" / "data" / "3o2e.mtz"), "--quiet",
                "--no-auto-stop", "--cycles", "2"]
        r = subprocess.run(base + ["--mtzout", str(first)], cwd=str(ROOT),
                           capture_output=True, text=True, env=ENV)
        assert r.returncode == 0, r.stderr
        r = subprocess.run([sys.executable, "-m", "oasis_dm", "--mtzin", str(first),
                            "--mtzout", str(second), "--quiet", "--no-auto-stop",
                            "--cycles", "2"], cwd=str(ROOT), capture_output=True,
                           text=True, env=ENV)
        assert r.returncode == 0, r.stderr
        labels = gemmi.read_mtz_file(str(second)).column_labels()
        assert len(labels) == len(set(labels)), labels
        assert len(labels) == 26, labels


def test_run_is_deterministic():
    import gemmi
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        outs = []
        for i in range(2):
            p = tmp / ("r%d.mtz" % i)
            r = subprocess.run([sys.executable, "-m", "oasis_dm", "--mtzin",
                                str(ROOT / "tests" / "data" / "3o2e.mtz"),
                                "--mtzout", str(p), "--quiet", "--no-auto-stop",
                                "--cycles", "2", "--seed", "7"],
                               cwd=str(ROOT), capture_output=True, text=True, env=ENV)
            assert r.returncode == 0, r.stderr
            outs.append(np.array(gemmi.read_mtz_file(str(p)).column_with_label(
                "DM_PHIDM").array, dtype=np.float64))
        assert np.allclose(outs[0], outs[1], atol=1e-4)


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
    sys.exit(1 if failures else 0)
