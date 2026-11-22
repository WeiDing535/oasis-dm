"""Phase-extension and anisotropic-correction tests (v2 engine)."""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import gemmi

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "tests" / "data"
ENV = dict(os.environ, PYTHONPATH=str(ROOT / "src"))


def _run(mtzin, mtzout, *extra):
    r = subprocess.run([sys.executable, "-m", "oasis_dm", "--mtzin", str(mtzin),
                        "--mtzout", str(mtzout), "--quiet"] + list(extra),
                       cwd=str(ROOT), capture_output=True, text=True, env=ENV)
    assert r.returncode == 0, r.stderr
    return Path(mtzout)


def _truth_rows(path):
    t = gemmi.read_mtz_file(str(DATA / "7e1d_ext.mtz"))
    labt = {c.label: i for i, c in enumerate(t.columns)}
    at = np.asarray(t.array)
    d = np.asarray(t.make_d_array(), dtype=np.float64)
    ext = d < 3.4
    idx = {tuple(r[:3]): i for i, r in enumerate(at)}
    return at, labt, ext, idx


def _score(path):
    at, labt, ext, idx = _truth_rows(path)
    m = gemmi.read_mtz_file(str(path))
    lab = {c.label: i for i, c in enumerate(m.columns)}
    a = np.asarray(m.array)
    rows = np.array([idx[tuple(r[:3])] for r in a])
    dphi = np.radians(at[rows, labt["PHIC"]] - a[:, lab["DM_PHIDM"]])
    fom = a[:, lab["DM_FOMDM"]]
    sel = ext[rows]
    cos = np.cos(dphi[sel])
    return cos.mean(), fom[sel].mean()


def test_extension_beats_no_extend():
    """Activating the unphased amplitudes must improve their true phases.

    Both runs pin the cycle count: with auto-stop the chosen cycle varies
    with the platform's FFT numerics and so does the measured gap.  The
    pinned gap is +0.012 cos on macOS, Linux and Windows alike.
    """
    with tempfile.TemporaryDirectory() as tmp:
        fixed = ["--no-auto-stop", "--cycles", "8"]
        on = _run(DATA / "7e1d_ext.mtz", Path(tmp) / "on.mtz", *fixed)
        off = _run(DATA / "7e1d_ext.mtz", Path(tmp) / "off.mtz",
                   "--no-extend", *fixed)
        cos_on, fom_on = _score(on)
        cos_off, fom_off = _score(off)
        # phases on the stripped reflections genuinely improve
        assert cos_on > cos_off + 0.01, (cos_on, cos_off)
        # the previously unphased reflections get a usable FOM
        assert fom_on > 0.2, fom_on
        # with extension off they are reported as unknown, as before
        assert fom_off == 0.0


def test_extension_output_columns_cover_all_reflections():
    with tempfile.TemporaryDirectory() as tmp:
        out = _run(DATA / "7e1d_ext.mtz", Path(tmp) / "o.mtz")
        m = gemmi.read_mtz_file(str(out))
        lab = {c.label: i for i, c in enumerate(m.columns)}
        a = np.asarray(m.array)
        src = gemmi.read_mtz_file(str(DATA / "7e1d_ext.mtz"))
        labs = {c.label: i for i, c in enumerate(src.columns)}
        s = np.asarray(src.array)
        n = len(a)
        assert n == len(s), (n, len(s))
        idx = {tuple(r[:3]): i for i, r in enumerate(s)}
        rows = np.array([idx[tuple(r[:3])] for r in a])
        # every reflection with a measured Fo now carries a phase and FOM
        fo = s[rows, labs["FOBS"]]
        has_fo = np.isfinite(fo) & (fo > 0)
        assert np.all(np.isfinite(a[has_fo, lab["DM_PHIDM"]]))
        assert np.all(a[has_fo, lab["DM_FOMDM"]] > 0.0)
        # FDM holds the measured amplitudes for those reflections
        assert np.allclose(a[has_fo, lab["DM_FDM"]], fo[has_fo], rtol=1e-3)


def test_extension_case_regenerates_with_tool():
    """tests/data/7e1d_ext.mtz must be reproducible via tools/make_ext_case.py."""
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "regen.mtz"
        r = subprocess.run([sys.executable, str(ROOT / "tools" / "make_ext_case.py"),
                            str(DATA / "7e1d.mtz"), str(out), "--dmin-phase", "3.4"],
                           cwd=str(ROOT), capture_output=True, text=True, env=ENV)
        assert r.returncode == 0, r.stderr
        a = np.asarray(gemmi.read_mtz_file(str(DATA / "7e1d_ext.mtz")).array)
        b = np.asarray(gemmi.read_mtz_file(str(out)).array)
        assert a.shape == b.shape
        assert np.allclose(np.nan_to_num(a), np.nan_to_num(b), equal_nan=True)


def test_aniso_is_safe_on_isotropic_data():
    """--aniso on the (isotropic) bundled case must not blow up the result."""
    with tempfile.TemporaryDirectory() as tmp:
        plain = _run(DATA / "7e1d.mtz", Path(tmp) / "plain.mtz")
        aniso = _run(DATA / "7e1d.mtz", Path(tmp) / "aniso.mtz", "--aniso")
        t = gemmi.read_mtz_file(str(DATA / "7e1d.mtz"))
        labt = {c.label: i for i, c in enumerate(t.columns)}
        at = np.asarray(t.array)
        idx = {tuple(r[:3]): i for i, r in enumerate(at)}
        for f in (plain, aniso):
            m = gemmi.read_mtz_file(str(f))
            lab = {c.label: i for i, c in enumerate(m.columns)}
            a = np.asarray(m.array)
            rows = np.array([idx[tuple(r[:3])] for r in a])
            dphi = np.radians(at[rows, labt["PHIC"]] - a[:, lab["DM_PHIDM"]])
            fom = a[:, lab["DM_FOMDM"]]
            w = fom
            yield_value = (w * np.cos(dphi)).sum() / w.sum()
            if f is plain:
                plain_wcos = yield_value
            else:
                assert yield_value > plain_wcos - 0.03, (yield_value, plain_wcos)


def test_aniso_recovers_synthetic_anisotropy():
    """--aniso must recover the damage a synthetic B along c* does."""
    with tempfile.TemporaryDirectory() as tmp:
        plain_iso = _run(DATA / "7e1d.mtz", Path(tmp) / "iso.mtz")   # ceiling
        off = _run(DATA / "7e1d_aniso.mtz", Path(tmp) / "off.mtz")
        on = _run(DATA / "7e1d_aniso.mtz", Path(tmp) / "on.mtz", "--aniso")

        def wcos(f, truth=DATA / "7e1d_aniso.mtz"):
            t = gemmi.read_mtz_file(str(truth))
            labt = {c.label: i for i, c in enumerate(t.columns)}
            at = np.asarray(t.array)
            idx = {tuple(r[:3]): i for i, r in enumerate(at)}
            m = gemmi.read_mtz_file(str(f))
            lab = {c.label: i for i, c in enumerate(m.columns)}
            a = np.asarray(m.array)
            rows = np.array([idx[tuple(r[:3])] for r in a])
            dphi = np.radians(at[rows, labt["PHIC"]] - a[:, lab["DM_PHIDM"]])
            fom = a[:, lab["DM_FOMDM"]]
            return float(np.cos(dphi).mean())

        iso, off_c, on_c = wcos(plain_iso), wcos(off), wcos(on)
        # the anisotropy genuinely hurts ...
        assert off_c < iso - 0.01, (off_c, iso)
        # ... and the correction recovers most of it
        assert on_c > off_c + (iso - off_c) * 0.6, (on_c, off_c, iso)


def test_aniso_case_regenerates_with_tool():
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "regen.mtz"
        r = subprocess.run([sys.executable, str(ROOT / "tools" / "make_aniso_case.py"),
                            str(DATA / "7e1d.mtz"), str(out)],
                           cwd=str(ROOT), capture_output=True, text=True, env=ENV)
        assert r.returncode == 0, r.stderr
        a = np.asarray(gemmi.read_mtz_file(str(DATA / "7e1d_aniso.mtz")).array)
        b = np.asarray(gemmi.read_mtz_file(str(out)).array)
        assert a.shape == b.shape
        assert np.allclose(np.nan_to_num(a), np.nan_to_num(b), equal_nan=True)


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


def test_difference_map_exports_complex_coefficients(tmp_path):
    output = _run(DATA / "6zpp.mtz", tmp_path / "difference.mtz",
                  "--map-update", "difference", "--no-auto-stop", "--cycles", "2")
    mtz = gemmi.read_mtz_file(str(output))

    def column(label):
        return np.asarray(mtz.column_with_label(label), dtype=float)

    def complex_coeff(amp, phase):
        return column(amp) * np.exp(1j * np.deg2rad(column(phase)))

    best = complex_coeff("DM_FWT", "DM_PHWT")
    difference = complex_coeff("DM_DELFWT", "DM_PHDELFWT")
    posterior = column("FOBS") * complex_coeff("DM_FOMDM", "DM_PHIDM")
    centric = mtz.spacegroup.operations().centric_flag_array(mtz.make_miller_array())
    measured = np.isfinite(column("FOBS")) & (column("FOBS") > 0)
    np.testing.assert_allclose(best[measured & centric], posterior[measured & centric],
                               rtol=1e-5, atol=1e-4)
    np.testing.assert_allclose(best[measured & ~centric] - difference[measured & ~centric],
                               posterior[measured & ~centric], rtol=1e-5, atol=1e-4)


def test_posterior_map_export_uses_observed_amplitude_scale(tmp_path):
    output = _run(DATA / "6zpp.mtz", tmp_path / "posterior.mtz", "--aniso",
                  "--map-update", "difference", "--map-output", "posterior",
                  "--no-auto-stop", "--cycles", "2")
    mtz = gemmi.read_mtz_file(str(output))
    fo = mtz.column_with_label("FOBS").array
    measured = np.isfinite(fo) & (fo > 0)
    fom = mtz.column_with_label("DM_FOMDM").array
    np.testing.assert_allclose(mtz.column_with_label("DM_FWT").array[measured],
                               (fo * fom)[measured], rtol=1e-6)
    np.testing.assert_allclose(mtz.column_with_label("DM_PHWT").array[measured],
                               mtz.column_with_label("DM_PHIDM").array[measured])
