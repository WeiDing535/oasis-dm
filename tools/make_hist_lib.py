#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate a protein-density histogram library in CCP4 ``hist.lib`` format.

``oasis_dm`` needs a resolution-dependent target distribution of protein electron
density for histogram matching.  CCP4 ships one as ``$CLIBD/hist.lib`` but that
file is not redistributable, so this script derives an independent library from
public PDB entries.

Method (per entry)
------------------
1. reduce the entry to amino-acid residues (waters, ligands, hydrogens dropped);
2. build the reflection list out to the finest target resolution and compute
   structure factors with ``gemmi.StructureFactorCalculatorX`` -- this carries the
   real atomic scattering factors f(s), which is what gives a density map its
   resolution-dependent contrast (a Gaussian atom model does not, and produces a
   target that is *worse* than no matching at all);
3. for each target resolution, truncate the reflection list and FFT onto a fixed
   grid, so the only thing that changes between columns is the truncation;
4. accumulate the histogram of the *protein* volume (constant-radius envelope)
   only.  Deposited entries differ wildly in solvent content -- and some ASUs
   contain only part of the biological assembly -- so averaging whole-cell
   histograms mixes unrelated solvent fractions.  The whole-cell target is
   assembled afterwards from the average protein distribution plus a canonical
   solvent fraction;
5. average the normalised histograms over entries and write 18 resolutions x 100
   bins in the same layout as ``hist.lib``.

Only the distribution *shape* is used by the matcher (the modified map is
re-standardised onto the current protein statistics), so the arbitrary overall
scale of the computed densities does not matter.

Usage
-----
    python tools/make_hist_lib.py --out src/oasis_dm/data/protein_hist.lib
    python tools/make_hist_lib.py --cache /tmp/pdb-cache      # keep downloads

PDB files are fetched from https://files.rcsb.org and cached locally.
"""
from __future__ import annotations

import argparse
import os
import sys
import urllib.request

import numpy as np
import gemmi

RESOLUTIONS = [10.0, 8.0, 7.0, 6.0, 5.0, 4.5, 4.0, 3.6, 3.3, 3.0,
               2.8, 2.6, 2.4, 2.2, 2.0, 1.8, 1.6, 1.5]
NBIN = 100
RMIN, RMAX = -1.0, 4.0
FINEST = RESOLUTIONS[-1]
SAMPLE = 2.5                # grid points per resolution element
MAX_POINTS = 3.0e6          # coarsen the grid beyond this to stay in memory
SOLC_REF = 0.5              # canonical solvent fraction of the assembled target
SOLV_SIGMA = 0.05           # solvent ripple width in protein-mean units
PROTEIN_RADIUS = 2.5        # A; constant-radius envelope around atoms

# Small, well-ordered, single-protein entries chosen for diversity of size and
# solvent content.  All are public PDB depositions.
STRUCTURES = ["1crn", "1ubq", "2vb1", "1ake"]


def fetch(pdb_id, cache):
    path = os.path.join(cache, pdb_id.upper() + ".cif")
    if not os.path.isfile(path):
        url = "https://files.rcsb.org/download/%s.cif" % pdb_id.upper()
        print("  downloading %s" % url)
        data = urllib.request.urlopen(url, timeout=60).read()
        with open(path, "wb") as fh:
            fh.write(data)
    return path


def protein_model(st):
    """First model, amino-acid residues only, hydrogens dropped."""
    model = st[0]
    model.remove_ligands_and_waters()
    model.remove_hydrogens()
    if not len(list(model.all())):
        raise ValueError('no amino-acid residues')
    return model


def hkl_sphere(cell, d_min):
    """All integer hkl with d >= d_min (P1 list, Friedel pairs included)."""
    rec = cell.reciprocal()
    # rows of rec.orth.mat are the reciprocal basis vectors a*, b*, c*
    A = np.array(rec.orth.mat, dtype=np.float64).reshape(3, 3)
    G = A @ A.T
    hm = int(np.ceil(cell.a / d_min)) + 1
    km = int(np.ceil(cell.b / d_min)) + 1
    lm = int(np.ceil(cell.c / d_min)) + 1
    h, k, l = np.meshgrid(np.arange(-hm, hm + 1), np.arange(-km, km + 1),
                          np.arange(-lm, lm + 1), indexing='ij')
    h, k, l = h.ravel(), k.ravel(), l.ravel()
    s2 = (G[0, 0] * h * h + G[1, 1] * k * k + G[2, 2] * l * l
          + 2 * G[0, 1] * h * k + 2 * G[0, 2] * h * l + 2 * G[1, 2] * k * l)
    # drop 000: without it the map oscillates about zero, so the protein region
    # has a positive mean and a finite spread -- the same convention as hist.lib,
    # whose axis runs from -1 to 4 (negative density inside the protein)
    keep = (s2 <= 1.0 / (d_min * d_min)) & (s2 > 0.0)
    return np.column_stack([h[keep], k[keep], l[keep]]).astype(np.int32), s2[keep]


def grid_size(cell, d_min, sample):
    return (int(np.ceil(cell.a / (d_min / sample))),
            int(np.ceil(cell.b / (d_min / sample))),
            int(np.ceil(cell.c / (d_min / sample))))


def density_from_sf(cell, hkl, sf, size):
    """rho = ifftn(conj(F)) * N/V -- the same convention oasis_dm uses."""
    nu, nv, nw = size
    F = np.zeros((nu, nv, nw), dtype=np.complex64)
    F[hkl[:, 0] % nu, hkl[:, 1] % nv, hkl[:, 2] % nw] = sf
    rho = np.fft.ifftn(F.conj()).real
    del F
    return rho * ((nu * nv * nw) / cell.volume)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', default='src/oasis_dm/data/protein_hist.lib')
    ap.add_argument('--cache', default=os.path.join(os.path.expanduser('~'),
                                                    '.cache', 'oasis-dm-pdb'))
    args = ap.parse_args(argv)
    os.makedirs(args.cache, exist_ok=True)

    acc = np.zeros((NBIN, len(RESOLUTIONS)), dtype=np.float64)
    used = []
    for pdb_id in STRUCTURES:
        try:
            st = gemmi.read_structure(fetch(pdb_id, args.cache))
            st.setup_entities()
            model = protein_model(st)
            cell = st.cell
            sg = gemmi.find_spacegroup_by_name(st.spacegroup_hm)
            sample = SAMPLE
            size = grid_size(cell, FINEST, sample)
            while float(np.prod(size)) > MAX_POINTS:
                sample *= 0.8
                size = grid_size(cell, FINEST, sample)
            print("  %s: %d atoms, grid %s" % (pdb_id, len(list(model.all())), size))

            env = gemmi.FloatGrid(size[0], size[1], size[2])
            env.set_unit_cell(cell)
            env.spacegroup = sg
            env.fill(0.0)
            env.mask_points_in_constant_radius(model, PROTEIN_RADIUS, 1.0,
                                               ignore_hydrogen=True)
            prot = np.array(env, dtype=np.float64) > 0.5
            frac = float(prot.mean())

            hkl, s2 = hkl_sphere(cell, FINEST)
            calc = gemmi.StructureFactorCalculatorX(cell)
            # gemmi computes one reflection at a time; ~20 us each, so a few
            # hundred thousand reflections is a few seconds per entry
            sf = np.array([calc.calculate_sf_from_model(model, h.tolist())
                           for h in hkl], dtype=np.complex64)
            print("       %d reflections, protein volume %.1f%%"
                  % (len(hkl), 100.0 * frac))

            for j, d in enumerate(RESOLUTIONS):
                print("       %5.2f A" % d, flush=True)
                sel = s2 <= 1.0 / (d * d)
                rho = density_from_sf(cell, hkl[sel], sf[sel], size)
                # hist.lib convention: mean protein density 1, solvent 0
                rho = rho / max(float(rho[prot].mean()), 1e-6)
                hist, _ = np.histogram(rho[prot], bins=NBIN, range=(RMIN, RMAX))
                acc[:, j] += hist.astype(np.float64) / max(float(hist.sum()), 1.0)
            used.append(pdb_id)
        except Exception as exc:                       # noqa: BLE001
            print("  %s: skipped (%s)" % (pdb_id, exc))

    if not used or acc.sum() == 0:
        raise SystemExit('no usable structures')

    axis = RMIN + (np.arange(NBIN) + 0.5) * (RMAX - RMIN) / NBIN
    solv = np.exp(-0.5 * (axis / SOLV_SIGMA) ** 2)
    solv /= solv.sum()
    cols, rave, rstd = [], [], []
    for j in range(len(RESOLUTIONS)):
        p = acc[:, j]
        p = p / p.sum()
        cell_pdf = (1.0 - SOLC_REF) * solv + SOLC_REF * p
        cell_pdf /= cell_pdf.sum()
        cols.append(cell_pdf)
        mean = float((axis * cell_pdf).sum())
        rave.append(mean)
        rstd.append(float(np.sqrt(((axis - mean) ** 2 * cell_pdf).sum())))

    out = os.path.abspath(args.out)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, 'w') as fh:
        fh.write("HIST %3d %3d\n" % (len(RESOLUTIONS), NBIN))
        fh.write("RESL" + "".join("%7.2f" % d for d in RESOLUTIONS) + "\n")
        fh.write("RMIN" + "".join("%7.3f" % RMIN for _ in RESOLUTIONS) + "\n")
        fh.write("RMAX" + "".join("%7.3f" % RMAX for _ in RESOLUTIONS) + "\n")
        fh.write("RAVE" + "".join("%7.3f" % v for v in rave) + "\n")
        fh.write("RSTD" + "".join("%7.3f" % v for v in rstd) + "\n")
        fh.write("DATA\n")
        for i in range(NBIN):
            fh.write("".join("%7.0f" % (1.0e4 * cols[j][i])
                             for j in range(len(RESOLUTIONS))) + "\n")
    print("wrote %s  (%d resolutions x %d bins from %s)"
          % (out, len(RESOLUTIONS), NBIN, ",".join(used)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
