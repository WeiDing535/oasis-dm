#!/usr/bin/env python3
"""Build a simulated-SAD test case (and NCS operator) from a PDB/mmCIF entry.

Generates the same column layout as the bundled test cases:

    H K L FOBS SIGFOBS FreeR PHIB FOM PHIC HLA HLB HLC HLD

* Fcalc / PHIC come from the model via gemmi's StructureFactorCalculatorX.
* Fo = |Fcalc| * lognormal amplitude noise (5%).
* The starting phase PHIB = PHIC + von Mises-ish noise (gaussian dphi with
  sigma chosen so the weighted mean phase error lands near --target-wmpe);
  FOM tracks the per-reflection cos of that error, floored at 0.05.
* FreeR marks ~5% of reflections as the free set with the value 0.

With --ncs-out the Kabsch rotation+translation between the first two
chains' CA atoms is written as a fractional NCS operator file for
`oasis_dm --ncs`.

Usage:
    python tools/make_test_case.py 3i0t.cif tests/data/3i0t.mtz \
        --dmin 2.5 --seed 3 --ncs-out tests/data/3i0t_ncs.txt
"""
import argparse
import sys
import os

import numpy as np
import gemmi


def kabsch(P, Q):
    """Rigid transform mapping P onto Q: returns (R, t) with Q ~ P @ R.T + t."""
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1.0, 1.0, d])
    R = (V @ D @ Wt).T          # row-vector convention: x' = x @ R.T + t
    t = Q.mean(0) - P.mean(0) @ R.T
    return R, t


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('model', help='PDB/mmCIF coordinates file')
    p.add_argument('mtzout')
    p.add_argument('--dmin', type=float, default=2.5)
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--target-wmpe', type=float, default=56.0,
                   help='starting weighted mean phase error (deg)')
    p.add_argument('--free-frac', type=float, default=0.05)
    p.add_argument('--ncs-out', default='',
                   help='also write the chain1->chain2 CA Kabsch operator here')
    a = p.parse_args(argv)

    st = gemmi.read_structure(a.model)
    st.setup_entities()
    st.remove_hydrogens()
    st.remove_alternative_conformations()
    cell, sg = st.cell, st.spacegroup_hm
    calc = gemmi.StructureFactorCalculatorX(cell)

    hkl = np.asarray(gemmi.make_miller_array(cell, gemmi.SpaceGroup(sg),
                                             a.dmin, unique=True),
                     dtype=np.float64)
    if len(hkl) < 500:
        p.error(f'only {len(hkl)} reflections at dmin={a.dmin}')

    Fc = np.array([calc.calculate_sf_from_model(st[0], list(map(int, h)))
                   for h in hkl])
    phi_c = np.rad2deg(np.angle(Fc)) % 360.0
    rng = np.random.default_rng(a.seed)
    fo = np.abs(Fc) * np.exp(rng.normal(0.0, 0.05, len(Fc)))
    sig = np.maximum(0.05 * fo, 0.5)

    # starting phases: PHIC + noise; FOM is a per-reflection error
    # estimate scattered around the SAD-like mean (independent of the draw,
    # so the weighted mean phase error lands near --target-wmpe).
    # Centric reflections get a sign error along their constraint axis --
    # an off-axis phase would be projected by the engine and would not
    # correspond to any real experiment.
    sigma = np.sqrt(-2.0 * np.log(np.cos(np.deg2rad(a.target_wmpe))))
    dphi = rng.normal(0.0, sigma, len(Fc))
    fom = np.clip(0.55 + 0.25 * rng.normal(0.0, 1.0, len(Fc)), 0.05, 0.95)

    phib = (phi_c + np.rad2deg(dphi)) % 360.0
    sgops = list(gemmi.SpaceGroup(sg).operations())
    n_cent = 0
    for i, hkl_row in enumerate(hkl):
        hh = [int(x) for x in hkl_row]
        for op in sgops:
            if list(op.apply_to_hkl(hh)) == [-hh[0], -hh[1], -hh[2]]:
                # sign error: a wrong-branch draw flips the phase by 180,
                # and the estimate weight tracks |cos(draw)|
                if np.cos(dphi[i]) < 0.0:
                    phib[i] = (phib[i] + 180.0) % 360.0
                fom[i] = np.clip(abs(np.cos(dphi[i])), 0.05, 0.95)
                dphi[i] = 0.0 if np.cos(dphi[i]) >= 0.0 else np.pi
                n_cent += 1
                break
    if n_cent:
        print(f'  centric reflections: {n_cent} (sign-error model)')
    # von Mises HL of a unimodal estimate: kappa from fom exactly the way
    # core.x_to_fom / core.fom_to_x invert each other (I1(k)/I0(k) = fom),
    # so the engine reads back the same figure of merit as the FOM column.
    from scipy.special import i0e, i1e
    from scipy.optimize import brentq

    def kappa_of(fom_val):
        return brentq(lambda k: i1e(k) / i0e(k) - fom_val, 1e-4, 200.0)

    kappa = np.array([kappa_of(f) for f in fom])
    hl = np.column_stack([kappa * np.cos(np.deg2rad(phib)),
                          kappa * np.sin(np.deg2rad(phib)),
                          np.zeros(len(Fc)), np.zeros(len(Fc))])

    free = np.where(rng.random(len(Fc)) < a.free_frac, 0, 19).astype(np.int32)

    out = gemmi.Mtz(with_base=True)
    out.spacegroup = gemmi.SpaceGroup(sg)
    out.cell = cell
    out.add_dataset('oasis')
    for lbl, typ in (('FOBS', 'F'), ('SIGFOBS', 'Q'), ('FreeR', 'I'),
                     ('PHIB', 'P'), ('FOM', 'W'), ('PHIC', 'P'),
                     ('HLA', 'A'), ('HLB', 'A'), ('HLC', 'A'), ('HLD', 'A')):
        out.add_column(lbl, typ)
    data = np.column_stack([hkl, fo, sig, free, phib, fom, phi_c, hl])
    out.set_data(np.ascontiguousarray(data))
    os.makedirs(os.path.dirname(os.path.abspath(a.mtzout)), exist_ok=True)
    out.write_to_file(a.mtzout)

    wmpe = np.rad2deg(np.arccos(np.clip(
        (fom * np.cos(dphi)).sum() / fom.sum(), -1, 1)))
    print(f'{a.mtzout}: {len(hkl)} reflections to {a.dmin} A, '
          f'{int((free == 0).sum())} free, starting wMPE {wmpe:.1f} deg')

    if a.ncs_out:
        chains = [ch for ch in st[0] if len(ch) > 20]
        if len(chains) < 2:
            p.error('--ncs-out needs at least two chains')
        ca = []
        for ch in chains:
            # residue seqid -> CA position: chains may start at different
            # residue numbers, so index-based pairing misaligns the Kabsch fit
            ca.append({r.seqid.num:
                       (r.find_atom('CA', '*').pos.x,
                        r.find_atom('CA', '*').pos.y,
                        r.find_atom('CA', '*').pos.z)
                       for r in ch if r.find_atom('CA', '*') is not None})
        if len(set(ca[0]) & set(ca[1])) < 20:
            p.error('fewer than 20 shared residues for the Kabsch fit')
        orth = np.asarray(cell.orth.mat.tolist())
        frac = np.asarray(cell.frac.mat.tolist())
        fvec = np.asarray(cell.frac.vec.tolist())
        with open(a.ncs_out, 'w') as fh:
            # operators chain0 -> every other chain: with conjugates under
            # the space group this averages the whole copy cluster
            for j in range(1, len(ca)):
                common = sorted(set(ca[0]) & set(ca[j]))
                if len(common) < 20:
                    continue
                P = np.array([ca[0][k] for k in common])
                Q = np.array([ca[j][k] for k in common])
                R, t = kabsch(P, Q)
                rmsd = float(np.sqrt(np.mean(np.sum(
                    (P @ R.T + t - Q) ** 2, axis=1))))
                ang = np.rad2deg(np.arccos(
                    np.clip((np.trace(R) - 1.0) / 2.0, -1, 1)))
                R_col = R.T                # kabsch is row-vector convention
                R_f = frac @ R_col @ orth  # column-vector fractional operator
                t_f = frac @ t + fvec
                # x' = R_f @ x + t_f  (column convention), file rows are R_f rows
                fh.write(f'# NCS operator {chains[0].name}->{chains[j].name} '
                         f'from {os.path.basename(a.model)} (Kabsch, CA)\n')
                fh.write(f'# rotation {ang:.1f} deg, CA rmsd {rmsd:.2f} A\n')
                for i in range(3):
                    fh.write(' '.join(f'{v:12.6f}' for v in
                                      (R_f[i, 0], R_f[i, 1], R_f[i, 2], t_f[i]))
                             + '\n')
                print(f'{a.ncs_out}: {chains[0].name}->{chains[j].name} '
                      f'rotation {ang:.1f} deg, CA rmsd {rmsd:.2f} A, '
                      f'{len(common)} CA pairs')
    return 0


if __name__ == '__main__':
    sys.exit(main())
