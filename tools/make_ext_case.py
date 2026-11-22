#!/usr/bin/env python3
"""Build a phase-extension test case from a phased MTZ.

Takes an MTZ with experimental phases (PHIB/FOM/HLA-HLD) and amplitudes
(FOBS) and blanks the phase columns for every reflection finer than
--dmin-phase, producing a file in which the coarse reflections carry the
starting phase information and the fine reflections have measured
amplitudes only.  Density modification must then extend the phases to the
full resolution of the amplitudes.

Usage:
    python3 tools/make_ext_case.py tests/data/7e1d.mtz tests/data/7e1d_ext.mtz \
        --dmin-phase 3.4

The truth phase column (default PHIC) and the free flag are left untouched
so the output can be scored with tools/evaluate.py.
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
    p.add_argument('--dmin-phase', type=float, default=3.4,
                   help='strip phases from reflections finer than this (A)')
    p.add_argument('--phib', default='PHIB')
    p.add_argument('--fom', default='FOM')
    p.add_argument('--hl', nargs=4, default=['HLA', 'HLB', 'HLC', 'HLD'])
    a = p.parse_args(argv)

    m = gemmi.read_mtz_file(a.mtzin)
    labels = [c.label for c in m.columns]
    for lbl in [a.phib, a.fom] + list(a.hl):
        if lbl not in labels:
            p.error(f'column {lbl} not found in {a.mtzin}')
    arr = np.array(m.array, copy=True)
    d = np.asarray(m.make_d_array(), dtype=np.float64)
    sel = d < a.dmin_phase
    for lbl in [a.phib, a.fom] + list(a.hl):
        arr[sel, labels.index(lbl)] = np.nan

    out = gemmi.Mtz(with_base=True)
    out.spacegroup = m.spacegroup
    out.cell = m.cell
    out.add_dataset('oasis_ext')
    for c in m.columns[3:]:
        out.add_column(c.label, c.type)
    out.set_data(np.ascontiguousarray(arr))
    os.makedirs(os.path.dirname(os.path.abspath(a.mtzout)), exist_ok=True)
    out.write_to_file(a.mtzout)
    n_strip = int(sel.sum())
    print(f'{a.mtzout}: {n_strip}/{len(d)} reflections stripped of phases '
          f'(d < {a.dmin_phase} A); {len(d) - n_strip} coarse reflections keep them')
    return 0


if __name__ == '__main__':
    sys.exit(main())
