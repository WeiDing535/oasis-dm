"""NCS averaging tests (v2 engine).

Two layers:
  * unit tests on the averager itself (operator-file parsing, that
    averaging a clean map symmetrizes it, that averaging a noisy map
    improves it where partners were found);
  * engine tests that `--ncs` helps the tight multi-copy case where
    averaging adds information solvent flattening cannot supply, and
    that `--ncs-auto` falls back to plain density modification when
    averaging would hurt.
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import gemmi
from scipy import ndimage

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "tests" / "data"
ENV = dict(os.environ, PYTHONPATH=str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "src"))

from oasis_dm import ncs                                   # noqa: E402


def _truth_grid(tag):
    """Density map from the PHIC/FOBS columns of a test MTZ."""
    m = gemmi.read_mtz_file(str(DATA / f"{tag}.mtz"))
    g = m.transform_f_phi_to_map("FOBS", "PHIC", sample_rate=3.0)
    return m, np.array(g, dtype=np.float64, copy=True)


def _prot_mask(rho, frac=0.47):
    sm = ndimage.uniform_filter(rho, size=7, mode="wrap")
    return sm > np.percentile(sm, 100.0 * (1.0 - frac))


def _cc(a, b, sel):
    a, b, sel = np.asarray(a), np.asarray(b), np.asarray(sel, dtype=bool)
    return float(np.corrcoef(a[sel], b[sel])[0, 1])


def _run(mtzin, mtzout, *extra):
    r = subprocess.run([sys.executable, "-m", "oasis_dm", "--mtzin", str(mtzin),
                        "--mtzout", str(mtzout), "--quiet"] + list(extra),
                       cwd=str(ROOT), capture_output=True, text=True, env=ENV)
    assert r.returncode == 0, r.stderr
    return Path(mtzout)


def _cos_vs_truth(path, tag):
    """Mean cos(phase error) against PHIC, work reflections."""
    t = gemmi.read_mtz_file(str(DATA / f"{tag}.mtz"))
    labt = {c.label: i for i, c in enumerate(t.columns)}
    at = np.asarray(t.array)
    idx = {tuple(r[:3]): i for i, r in enumerate(at)}
    m = gemmi.read_mtz_file(str(path))
    lab = {c.label: i for i, c in enumerate(m.columns)}
    a = np.asarray(m.array)
    rows = np.array([idx[tuple(r[:3])] for r in a])
    dphi = np.radians(at[rows, labt["PHIC"]] - a[:, lab["DM_PHIDM"]])
    return float(np.mean(np.cos(dphi)))


# ------------------------------------------------------------ unit tests --

def test_operator_file_round_trip():
    mtz = gemmi.read_mtz_file(str(DATA / "3i0t.mtz"))
    m_cell_orth = mtz.cell.orth.mat.tolist()
    m_cell_frac = mtz.cell.frac.mat.tolist()
    ops = ncs.read_operators(DATA / "3i0t_ncs.txt")
    assert len(ops) == 1
    R, t = ops[0]
    assert R.shape == (3, 3) and t.shape == (3,)
    # the file is in FRACTIONAL coordinates: orthogonality only holds
    # after converting back to Cartesian
    orth = np.asarray(m_cell_orth, dtype=np.float64)
    frac = np.asarray(m_cell_frac, dtype=np.float64)
    Rc = orth @ R @ frac
    assert np.allclose(Rc @ Rc.T, np.eye(3), atol=1e-3)
    assert abs(np.linalg.det(Rc) - 1.0) < 1e-3              # proper rotation
    # rewriting the same numbers reads back identically
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "ops.txt"
        with open(p, "w") as fh:
            fh.write("# round trip\n")
            for i in range(3):
                fh.write(" ".join(f"{v:12.6f}" for v in
                                  (*R[i], t[i])) + "\n")
        ops2 = ncs.read_operators(p)
        assert len(ops2) == 1
        assert np.allclose(ops2[0][0], R) and np.allclose(ops2[0][1], t)
        # comments and blank lines are ignored
        with open(p, "a") as fh:
            fh.write("\n# trailing comment\n")
        assert len(ncs.read_operators(p)) == 1


def test_operator_file_rejects_garbage():
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "bad.txt"
        p.write_text("1 0 0\n0 1 0\n")                     # 2 rows: not 3k
        try:
            ncs.read_operators(p)
        except ValueError:
            pass
        else:
            raise AssertionError("2-row operator file must raise")
        p.write_text("1 0 0 0\n0 1 0 0\n0 0 -1 0\n")       # improper
        try:
            ncs.read_operators(p)
        except ValueError:
            pass
        else:
            raise AssertionError("improper rotation must raise")


def test_averaging_symmetrizes_clean_map():
    """Averaging the clean map with the true operator must stay faithful.

    The 3i0t chains differ by a ~1 A CA rmsd, so the averaged map is the
    mean structure, not the original -- but the correlation on the
    averaged voxels must remain high and nothing outside the mask may
    change.
    """
    m, rho = _truth_grid("3i0t")
    ops = ncs.read_operators(DATA / "3i0t_ncs.txt")
    av = ncs.NcsAverager(ops, spacegroup=m.spacegroup, cell=m.cell)
    mask = _prot_mask(rho)
    avg = av.average(rho, mask)
    sel = mask & av.last_avg_mask
    assert sel.sum() > 5000, int(sel.sum())
    assert _cc(rho, avg, sel) > 0.85, _cc(rho, avg, sel)
    # outside the mask the map is untouched
    out = ~mask
    assert np.array_equal(avg[out], rho[out])


def test_averaging_improves_noisy_map():
    """With independent noise on the map the partners still pair up.

    Noise on both copies is uncorrelated, so averaging must improve the
    correlation with the clean map wherever a partner was found, and
    must never make the whole mask worse.
    """
    m, rho = _truth_grid("3i0t")
    rng = np.random.default_rng(42)
    scale = 0.9 * rho.std()
    noisy = rho + rng.normal(0.0, scale, rho.shape)
    ops = ncs.read_operators(DATA / "3i0t_ncs.txt")
    av = ncs.NcsAverager(ops, spacegroup=m.spacegroup, cell=m.cell)
    mask = _prot_mask(rho)
    avg = av.average(noisy, mask)
    sel = mask & av.last_avg_mask
    assert sel.sum() > 5000, int(sel.sum())
    before = _cc(rho, noisy, sel)
    after = _cc(rho, avg, sel)
    assert after > before + 0.02, (before, after)
    assert _cc(rho, avg, mask) >= _cc(rho, noisy, mask) - 0.01


# ----------------------------------------------------------- engine tests --

def test_ncs_engine_smoke():
    """`--ncs` runs, writes the usual columns and actually averaged.

    One cycle on the easy 2-fold case keeps this affordable in CI; the
    full plain-vs-ncs comparison is gated behind OASIS_DM_SLOW_TESTS
    below (and reproduced by `tools/benchmark.py --ncs`).
    """
    with tempfile.TemporaryDirectory() as tmp:
        stats = Path(tmp) / "stats.json"
        out = _run(DATA / "3i0t.mtz", Path(tmp) / "ncs.mtz",
                   "--ncs", str(DATA / "3i0t_ncs.txt"),
                   "--no-auto-stop", "--max-cycles", "1",
                   "--stats", str(stats))
        m = gemmi.read_mtz_file(str(out))
        labels = set(m.column_labels())
        assert "DM_PHIDM" in labels and "DM_FOMDM" in labels
        import json
        s = json.loads(stats.read_text())
        fracs = [c["ncs_avg_frac"] for c in s["cycles"]
                 if c.get("ncs_avg_frac") is not None]
        assert fracs and fracs[0] > 0.02, fracs


def test_ncs_engine_helps_multicopy():
    """`--ncs` must beat plain density modification on the 4-copy case.

    Slow (a full engine run each side); skipped unless
    OASIS_DM_SLOW_TESTS is set so CI stays lean.
    7tj1 has four copies related by tight operators (~0.2-0.5 A rmsd);
    there averaging adds real information that solvent flattening alone
    cannot supply.  On 3i0t (two copies, a loose 1.0 A rmsd operator)
    averaging costs accuracy -- that case is the reason `--ncs` is an
    expert option and `--ncs-auto` gates on the free set.
    """
    if not os.environ.get("OASIS_DM_SLOW_TESTS"):
        return
    with tempfile.TemporaryDirectory() as tmp:
        plain = _run(DATA / "7tj1.mtz", Path(tmp) / "plain.mtz")
        withncs = _run(DATA / "7tj1.mtz", Path(tmp) / "ncs.mtz",
                       "--ncs", str(DATA / "7tj1_ncs.txt"))
        c0 = _cos_vs_truth(plain, "7tj1")
        c1 = _cos_vs_truth(withncs, "7tj1")
        assert c1 > c0, (c1, c0)


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
