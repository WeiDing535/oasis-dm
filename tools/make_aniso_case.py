#!/usr/bin/env python3
"""Build an anisotropic-amplitude test case from a phased MTZ.

Applies a synthetic direction-dependent Debye-Waller factor to the
amplitudes (default: an extra B of 60 A^2 along c*, compensated along a*
and b*), keeping the phases untouched.  The result is a dataset whose
*phases* are as good as the input but whose *amplitudes* fall off
anisotropically -- the situation --aniso is designed to correct.

Usage:
    python3 tools/make_aniso_case.py tests/data/7e1d.mtz tests/data/7e1d_aniso.mtz

Only orthogonal cells are handled (the bundled test cases are).
"""
import argparse
import sys
import os

import numpy as np
import gemmi


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('mtzin')
    p.add_argument('mtzout')
    p.add_argument('--b-cstar', type=float, default=60.0,
                   help='extra B along c* (A^2); a*/b* get -B/2 (default 60)')
    p.add_argument('--fo', default='FOBS')
    p.add_argument('--sigfo', default='SIGFOBS')
    a = p.parse_args(argv)

    m = gemmi.read_mtz_file(a.mtzin)
    cell = m.cell
    if abs(cell.alpha - 90) > 1e-6 or abs(cell.beta - 90) > 1e-6 \
            or abs(cell.gamma - 90) > 1e-6:
        p.error('non-orthogonal cell: this generator only handles orthogonal cells')
    labels = [c.label for c in m.columns]
    for lbl in (a.fo, a.sigfo):
        if lbl not in labels:
            p.error(f'column {lbl} not found in {a.mtzin}')
    arr = np.asarray(m.array, copy=True)
    hkl = np.asarray(m.make_miller_array(), dtype=np.float64)
    qz = hkl[:, 2] / cell.c
    qxy2 = (hkl[:, 0] / cell.a) ** 2 + (hkl[:, 1] / cell.b) ** 2
    fac = np.exp(-a.b_cstar / 4.0 * qz ** 2 + a.b_cstar / 8.0 * qxy2)
    for lbl in (a.fo, a.sigfo):
        arr[:, labels.index(lbl)] *= fac

    out = gemmi.Mtz(with_base=True)
    out.spacegroup = m.spacegroup
    out.cell = m.cell
    out.add_dataset('oasis_aniso')
    for c in m.columns[3:]:
        out.add_column(c.label, c.type)
    out.set_data(np.ascontiguousarray(arr))
    os.makedirs(os.path.dirname(os.path.abspath(a.mtzout)), exist_ok=True)
    out.write_to_file(a.mtzout)
    print(f'{a.mtzout}: amplitudes scaled by {fac.min():.3f}-{fac.max():.3f} '
          f'(B_c* = +{a.b_cstar}, B_a*/B_b* = -{a.b_cstar / 2} A^2)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
