"""Non-crystallographic symmetry: operator I/O, map averaging, auto-detection.

Conventions
-----------
* An operator is a fractional-coordinate rigid motion ``x' = R @ x + t``
  (``R`` a proper rotation, det = +1).  Fractional coordinates are the
  natural system: they are what the density grid is indexed in.
* The density grid ``rho[i, j, k]`` holds the value at fractional position
  ``(i/N1, j/N2, k/N3)`` and is *periodic* (it is the inverse FFT of the
  full P1 reciprocal grid), so rotations/translations wrap at the cell
  edges exactly as the physics does -- no boundary artifacts.
* Averaging is orbit-based: a point is averaged with the images under the
  operator's forward powers (and, for open operators that never close on
  the identity, its inverse), wherever both the point and the image lie in
  the protein mask.

Auto-detection
--------------
A coarse SO(3) Euler scan (20 deg) of a downsampled protein mask, scored by
the translation-FFT cross-correlation peak between the mask and its rotated
copy.  Crystallographic rotations are filtered out before scoring (they
would trivially score 1.0).  The best candidate is refined at 5 deg, then
its translation is refined on the full-resolution mask with one more FFT
correlation.
"""
import itertools

import numpy as np
from scipy import ndimage

__all__ = ['read_operators', 'NcsAverager', 'detect_ncs', 'rot_zxz',
           'cart_to_frac']


# --------------------------------------------------------------- operators --

def cart_to_frac(R_c, cell):
    """Cartesian rotation -> fractional rotation (gemmi matrices)."""
    orth = np.asarray(cell.orth.mat.tolist(), dtype=np.float64)
    frac = np.asarray(cell.frac.mat.tolist(), dtype=np.float64)
    return frac @ R_c @ orth


def _idx_matrix(R_f, t_f, shape):
    """Fractional (R, t) -> scipy index-space (matrix, offset).

    scipy's affine_transform evaluates output[p] = input[matrix @ p +
    offset]; with M[i,j] = R_f[i,j] * N_i / N_j and offset = t_f * N this
    samples the density at the fractional image R_f @ x + t of grid point p.
    """
    n = np.asarray(shape, dtype=np.float64)
    M = R_f * (n[:, None] / n[None, :])
    return M, t_f * n


def _frac_metric(cell):
    """Metric tensor of fractional coordinates (cartesian inner product)."""
    orth = np.asarray(cell.orth.mat.tolist(), dtype=np.float64)
    return orth.T @ orth


def _is_proper(R_f, cell=None, tol=2e-3):
    """Proper rotation in *fractional* coordinates.

    A fractional rotation is not orthonormal under the plain inner product;
    it preserves the lattice metric G = orth^T orth.  Without a cell only
    the determinant can be checked.
    """
    if abs(np.linalg.det(R_f) - 1.0) > 5e-3:
        return False
    if cell is None:
        return True
    G = _frac_metric(cell)
    return np.allclose(R_f.T @ G @ R_f, G, atol=tol * np.trace(G))


def read_operators(path, cell=None):
    """Read NCS operators from a text file.

    Each operator is three lines of four numbers (a row of R followed by
    the translation component)::

        r11 r12 r13 t1
        r21 r22 r23 t2
        r31 r32 r33 t3

    Blank lines and #-comments are ignored; 12 numbers on one line also
    work.  Returns a list of (R_f, t_f) tuples.
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.split('#')[0].strip()
            if line:
                rows.append([float(x) for x in line.split()])
    if rows and all(len(r) == 12 for r in rows):
        # single-line-per-operator form: unpack into r-r-r-t rows
        groups = []
        for r in rows:
            for i in range(3):
                groups.append(r[4 * i:4 * i + 3] + [r[4 * i + 3]])
    elif rows and all(len(r) == 4 for r in rows):
        groups = rows
    else:
        raise ValueError(
            f'{path}: expected 3 lines of "r r r t" per operator '
            f'(or 12 numbers on one line); got '
            f'{[len(r) for r in rows][:6]}...')
    if len(groups) == 0 or len(groups) % 3 != 0:
        raise ValueError(f'{path}: expected 3 rows per operator, got '
                         f'{len(groups)}')
    ops = []
    for k in range(len(groups) // 3):
        g = groups[3 * k:3 * k + 3]
        M = np.asarray([r[:3] for r in g], dtype=np.float64)
        t = np.asarray([r[3] for r in g], dtype=np.float64)
        if not _is_proper(M, cell):
            raise ValueError(
                f'{path}: operator {k + 1} is not a proper rotation '
                f'(det={np.linalg.det(M):.4f}); improper NCS is not supported')
        ops.append((M, t))
    return ops


# --------------------------------------------------------------- averaging --

class NcsAverager:
    """Applies orbit averaging of a density map inside the NCS mask.

    When a space group is supplied the averaging *modification* is
    propagated to every crystallographic copy of the asymmetric unit.
    The orbit averaging itself is only valid on the tile of the cell
    where the operator maps protein onto protein (the mask intersection
    selects it); without the propagation the modified map would break
    crystallographic symmetry and the structure factors extracted from
    it would be corrupted by the asymmetric residue.
    """

    def __init__(self, ops, max_order=6, spacegroup=None, cell=None,
                 cc_thresh=0.25):
        self.ops = list(ops)
        self.max_order = max_order
        self.spacegroup = spacegroup
        self.cell = cell
        self.cc_thresh = cc_thresh
        self.last_frac = None   # fraction of grid points actually averaged
        self._cc_win = None     # cached local-correlation window (voxels)

    def _close_to_identity(self, R, t):
        """Is the operator (R, t) the identity up to a lattice vector?

        A two-fold NCS operator estimated at 179.5 deg squares to a
        ~1 deg rotation plus a lattice shift; that orbit *does* close.
        The tolerance is therefore set in cartesian angstroms, not in
        matrix elements: a partner within ~2.5 A of the point itself is
        a degenerate self-average and must not enter the orbit.
        """
        if not np.allclose(R, np.eye(3), atol=0.05):    # ~3 deg
            return False
        tf = t - np.round(t)
        if self.cell is None:
            return np.allclose(tf, 0.0, atol=5e-2)
        orth = np.asarray(self.cell.orth.mat.tolist(), dtype=np.float64)
        return float(np.linalg.norm(tf @ orth.T)) < 2.5

    def _sg_ops(self):
        """Space-group operators as fractional (R, t) pairs."""
        if self.spacegroup is None:
            return []
        out = []
        for op in self.spacegroup.operations():
            R = np.asarray(op.rot, dtype=np.float64) / op.DEN
            t = np.asarray(op.tran, dtype=np.float64) / op.DEN
            out.append((R, t))
        return out

    def _same_op(self, a, b):
        """Do two operators agree up to a lattice translation?"""
        Ra, ta = a
        Rb, tb = b
        if not np.allclose(Ra, Rb, atol=1e-3):
            return False
        d = (ta - tb)
        d = d - np.round(d)
        if self.cell is None:
            return np.allclose(d, 0.0, atol=1e-2)
        orth = np.asarray(self.cell.orth.mat.tolist(), dtype=np.float64)
        return float(np.linalg.norm(d @ orth.T)) < 0.5

    def _op_set(self):
        """Operators closed under space-group conjugation.

        A non-crystallographic operator O is a true symmetry of the
        *content* only on the tile of the cell it was fitted to.  The
        crystallographic images of that tile are related by the
        conjugates s_j O s_j^-1, each valid on its own tile, and the
        whole set {s_j O s_j^-1} is closed under conjugation -- so
        averaging with the full set produces a map that is exactly
        space-group symmetric by construction.  Conjugates that are
        themselves crystallographic operators are dropped (they would
        average each point with an identical copy of itself).
        """
        ops = list(self.ops)
        sg = self._sg_ops()
        if not sg:
            return ops
        out = []
        for R, t in ops:
            for S, u in sg:
                Si = np.linalg.inv(S)
                Rc = S @ R @ Si
                tc = -Rc @ u + S @ t + u
                cand = (Rc, tc)
                if self._close_to_identity(*cand):
                    continue
                if any(self._same_op(cand, x) for x in sg):
                    continue        # crystallographic: nothing to average
                if any(self._same_op(cand, x) for x in out):
                    continue        # duplicate up to a lattice vector
                out.append(cand)
        return out

    def _local_cc_window(self, shape):
        """Window (voxels) for the local-correlation test, ~8 A wide."""
        if self._cc_win is not None and self._cc_win[0].shape == shape:
            return self._cc_win[1]
        if self.cell is None:
            win = (5, 5, 5)
        else:
            orth = np.asarray(self.cell.orth.mat.tolist(), dtype=np.float64)
            win = []
            for i in range(3):
                axis_len = float(np.linalg.norm(orth[:, i]))
                win.append(int(round(8.0 / (axis_len / shape[i]))))
            win = tuple(max(3, min(15, w)) for w in win)
        self._cc_win = (np.empty(shape, dtype=bool), win)
        return win

    def _local_cc(self, rho, rot, rstats=None):
        """Local correlation between a map and its rotated image.

        The plain mask intersection is a weak partner filter when the
        protein fills about half the cell (random overlap of the mask
        with itself is then comparable to the true NCS overlap).  The
        local correlation distinguishes genuine NCS mates (positive)
        from accidental protein-protein overlaps (near zero).

        ``rstats`` optionally carries ``(uf(rho), uf(rho*rho))`` cached
        by :meth:`average` -- the map side is identical for every
        operator variant and need not be recomputed.
        """
        win = self._local_cc_window(rho.shape)
        uf = lambda a: ndimage.uniform_filter(a, size=win, mode='wrap')
        if rstats is None:
            rstats = (uf(rho), uf(rho * rho))
        mr, mr2 = rstats
        vr = mr2 - mr * mr
        mo = uf(rot)
        vo = uf(rot * rot) - mo * mo
        cov = uf(rho * rot) - mr * mo
        return cov / np.sqrt(np.maximum(vr * vo, 1e-12))

    def _powers(self, R, t):
        """Forward powers of the operator until it closes on the identity.

        Returns (powers, closed): powers is a list of (R_k, t_k) with
        Op^k = R^k, t_k = (I + R + ... + R^{k-1}) t.  Closure is judged
        in cartesian angstroms (see ``_close_to_identity``): a two-fold
        estimated at 179.5 deg closes at order 2, its residual ~1 deg
        being operator-estimation noise rather than a real third copy.
        """
        powers = []
        Rk, tk = R, t
        for _ in range(self.max_order):
            if self._close_to_identity(Rk, tk):
                return powers, True      # orbit closed before this power
            powers.append((Rk, tk))
            Rk = R @ Rk
            tk = R @ tk + t
        return powers[:1], False         # open operator: pairwise only

    def _lattice_variants(self):
        """Fractional lattice shifts to try for one operator.

        An operator fitted to unwrapped fractional coordinates obeys
        ``O(f + L) = O(f) + R @ L``.  On the periodic grid only wrapped
        coordinates are available, and because ``R`` is not an integer
        matrix, applying ``O`` to the wrapped representative of a point
        is NOT equivalent to applying it to the representative the fit
        used.  Molecules that straddle a unit-cell boundary in
        fractional space (very common) would silently receive partners
        from the wrong image.  Trying every small lattice shift and
        keeping, per voxel, the locally best-correlated partner removes
        that ambiguity.
        """
        return [np.asarray(L, dtype=np.float64)
                for L in itertools.product((-1, 0, 1), repeat=3)]

    def _accumulate(self, rho, mask, R, t, out, wsum, rstats=None):
        shape = rho.shape
        powers, closed = self._powers(R, t)
        seq = list(powers)
        if not closed:
            # open operator: also walk backwards so both copies average,
            # skipping degenerate near-self steps
            Ri, ti = R.T, -R.T @ t
            for _ in range(self.max_order):
                if self._close_to_identity(Ri, ti):
                    break
                seq.append((Ri, ti))
                Ri = R.T @ Ri
                ti = R.T @ ti - R.T @ t
        maskf = mask.astype(np.float64)
        min_ov = max(30, 0.002 * float(mask.sum()))
        win = self._local_cc_window(shape)
        for Rk, tk in seq:
            # pass 1 (cheap): rotate the mask under every lattice variant
            # and keep only variants whose domain overlaps the mask.
            # (The variants cannot be reduced to grid rolls of one warp:
            # the roll identity rot_L(x) = rot_0(x + n*L) breaks where
            # x + n*L leaves the box, because rho is periodic on the
            # integer lattice, not on R @ (integer lattice).)
            cand = []
            for L in self._lattice_variants():
                tL = tk + Rk @ L
                M, off = _idx_matrix(Rk, tL, shape)
                mrot = ndimage.affine_transform(
                    maskf, M, offset=off, order=1,
                    mode='grid-wrap') > 0.35
                if int(np.count_nonzero(mask & mrot)) >= min_ov:
                    cand.append((M, off, mrot))
            if not cand:
                continue
            # pass 2: local correlation of the surviving variants
            stats = []
            for M, off, mrot in cand:
                rot = ndimage.affine_transform(
                    rho, M, offset=off, order=1, mode='grid-wrap')
                lcc = self._local_cc(rho, rot, rstats)
                stats.append((mrot, rot, lcc))
            # pass 3: assign each voxel the lattice variant of its nearest
            # high-confidence seed voxel, then blend the partner values in
            # with a soft weight.
            #
            # The protein mask is often ONE connected component (the NCS
            # copies touch each other or are merged by the mask blur), so
            # a per-component vote cannot separate the copies: every copy
            # needs a different lattice variant, and a whole-mask vote
            # would pick a compromise variant that is right for no one.
            # Instead: voxels whose smoothed local cc clearly favours one
            # variant (above ``0.30``, far beyond the noise of a second
            # box-smoothing) act as seeds, and every other mask voxel
            # adopts the seed variant of its Voronoi cell -- a partition
            # of the mask that follows the molecular copies.  The blend
            # weight is the double-smoothed local cc of the *assigned*
            # variant (not the max over variants, which noise would
            # inflate): lcc_smooth <= 0.10 -> no blend,
            # lcc_smooth >= 0.25 -> full replacement.
            LCC = np.stack([s[2] for s in stats]).astype(np.float32)
            D = np.empty_like(LCC)
            for v in range(LCC.shape[0]):
                D[v] = ndimage.uniform_filter(
                    LCC[v], size=win, mode='wrap')
            best_idx = np.argmax(D, axis=0)
            Dbest = np.max(D, axis=0)
            mrot_best = np.zeros(shape, dtype=bool)
            for v in range(LCC.shape[0]):
                mrot_best |= (best_idx == v) & stats[v][0]
            seeds = mask & (Dbest > 0.30) & mrot_best
            if int(np.count_nonzero(seeds)) < 200:
                seeds = mask & (Dbest > 0.20) & mrot_best
            if int(np.count_nonzero(seeds)) < 200:
                continue        # no confident partner anywhere: skip op
            _, (ii, jj, kk) = ndimage.distance_transform_edt(
                ~seeds, return_indices=True)
            vidx = best_idx[ii, jj, kk]
            # weight = the DOUBLE-smoothed local cc of the assigned
            # variant: within a Voronoi cell whose seed chose it with
            # Dbest > 0.30 the whole cell averages at full weight (the
            # symmetrization gain is proportional to averaging strength),
            # while cells whose correlation never rose above the noise
            # stay near zero.  Unlike a max over variants the smoothed
            # assigned value is not inflated by winner selection.
            for v in range(LCC.shape[0]):
                sel = mask & (vidx == v) & stats[v][0]
                if not sel.any():
                    continue
                w = np.clip((D[v] - 0.10) / 0.15, 0.0, 1.0) * sel
                wsum = wsum + w
                out = out + w * stats[v][1]
        return out, wsum

    def average(self, rho, mask):
        """Average ``rho`` over the NCS orbits inside boolean ``mask``.

        The operator set is closed under space-group conjugation, so
        every crystallographic copy of the NCS-related pair is averaged
        with its own mate and the result stays space-group symmetric.
        """
        rho = np.asarray(rho, dtype=np.float64)
        win = self._local_cc_window(rho.shape)
        uf = lambda a: ndimage.uniform_filter(a, size=win, mode='wrap')
        rstats = (uf(rho), uf(rho * rho))
        out = rho.copy()
        wsum = np.ones(rho.shape, dtype=np.float64)
        for R, t in self._op_set():
            out, wsum = self._accumulate(rho, mask, R, t, out, wsum, rstats)
        self.last_frac = float(np.mean(wsum > 1.5))
        self.last_avg_mask = wsum > 1.5      # voxels that gained a partner
        return np.where(wsum > 1.0, out / np.maximum(wsum, 1.0), rho)

    def wrap(self, modify, mask):
        """Compose NCS averaging after a density-modification function."""
        def modify_with_ncs(x):
            return self.average(modify(x), mask)
        return modify_with_ncs


# ----------------------------------------------------------- auto-detection --

def rot_zxz(alpha, beta, gamma):
    """Euler angles (rad, Z-X-Z) -> Cartesian rotation matrix."""
    ca, sa = np.cos(alpha), np.sin(alpha)
    cb, sb = np.cos(beta), np.sin(beta)
    cg, sg = np.cos(gamma), np.sin(gamma)
    Rz = lambda c, s: np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])
    Rx = lambda c, s: np.array([[1.0, 0, 0], [0, c, -s], [0, s, c]])
    return Rz(ca, sa) @ Rx(cb, sb) @ Rz(cg, sg)


def _sg_rotations(spacegroup):
    """Distinct rotation parts of the space-group operations."""
    rots = [np.eye(3)]
    for op in spacegroup.operations():
        r = np.asarray(op.rot, dtype=np.float64) / op.DEN
        # NB: np.any(<generator>) is always True -- use the builtin any()
        if not any(np.allclose(r, x, atol=1e-6) for x in rots):
            rots.append(r)
    return rots


def _fft_corr_peak(a, b):
    """Peak of the circular cross-correlation sum_p a[p] b[p+t]."""
    A = np.fft.rfftn(a)
    B = np.fft.rfftn(b)
    c = np.fft.irfftn(np.conj(A) * B, s=a.shape)
    t = np.unravel_index(np.argmax(c), c.shape)
    return float(c[t]), t


def _voxel_size(cell, shape):
    """Cell axis lengths / grid dimensions (A per voxel, per axis)."""
    orth = np.asarray(cell.orth.mat.tolist(), dtype=np.float64)
    return np.array([np.linalg.norm(orth[:, i]) / shape[i]
                     for i in range(3)])


def _rot_angle_between(A, B):
    """Geodesic angle (rad) between two rotation matrices."""
    c = (np.trace(np.asarray(A).T @ np.asarray(B)) - 1.0) / 2.0
    return float(np.arccos(np.clip(c, -1.0, 1.0)))


def _mean_local_cc(rho, core, R_f, t_f, win, mrot=None, variants=False):
    """Mean local correlation between ``rho`` and its image under (R_f, t_f),
    evaluated over ``core`` intersected with the image of ``core``.

    With ``variants=True`` the score is the per-voxel maximum over the 27
    lattice variants ``t + R @ L`` of the operator (see
    ``NcsAverager._lattice_variants``): molecules that straddle a
    unit-cell boundary in fractional coordinates have partners under a
    non-mod-1 representative, and a single wrapped operator would
    underestimate a genuinely good operator.  ``mrot`` is ignored in
    that mode.
    """
    shape = rho.shape
    uf = lambda q: ndimage.uniform_filter(q, size=win, mode='wrap')
    if not variants:
        M, off = _idx_matrix(R_f, t_f, shape)
        rot = ndimage.affine_transform(
            rho, M, offset=off, order=1, mode='grid-wrap')
        if mrot is None:
            mrot = ndimage.affine_transform(
                core.astype(np.float64), M, offset=off, order=1,
                mode='grid-wrap') > 0.35
        reg = core & mrot
        if float(reg.sum()) < 30:
            return -1.0, reg
        m1, m2 = uf(rho), uf(rot)
        v1 = uf(rho * rho) - m1 * m1
        v2 = uf(rot * rot) - m2 * m2
        cov = uf(rho * rot) - m1 * m2
        lcc = cov / np.sqrt(np.maximum(v1 * v2, 1e-12))
        return float(np.mean(lcc[reg])), reg
    # variant-aware scoring
    m1 = uf(rho)
    v1 = uf(rho * rho) - m1 * m1
    coref = core.astype(np.float64)
    best = None
    reg_any = np.zeros(shape, dtype=bool)
    for L in itertools.product((-1, 0, 1), repeat=3):
        tL = t_f + R_f @ np.asarray(L, dtype=np.float64)
        M, off = _idx_matrix(R_f, tL, shape)
        mrot = ndimage.affine_transform(
            coref, M, offset=off, order=1, mode='grid-wrap') > 0.35
        sel = core & mrot
        if float(sel.sum()) < 30:
            continue
        rot = ndimage.affine_transform(
            rho, M, offset=off, order=1, mode='grid-wrap')
        m2 = uf(rot)
        v2 = uf(rot * rot) - m2 * m2
        cov = uf(rho * rot) - m1 * m2
        lcc = cov / np.sqrt(np.maximum(v1 * v2, 1e-12))
        cand = np.where(sel, lcc, -1.0)
        if best is None:
            best = cand
        else:
            best = np.maximum(best, cand)
        reg_any |= sel
    if best is None or float(reg_any.sum()) < 30:
        return -1.0, reg_any
    return float(np.mean(best[reg_any])), reg_any


def _nm_refine(rho_d, core_d, cell, x0, sg_cart, win, maxiter=350):
    """Nelder-Mead refinement of (rotation vector, translation).

    ``x0`` = [w (3, cartesian rotvec), t (3, fractional)].  The objective
    is the mean *local* correlation between the map and its image over
    the core region -- the same quantity the averager relies on.  A
    global Pearson score would be dominated by near-crystallographic
    operators (the map is space-group symmetric, so those correlate by
    construction), and self-correlation would drift to the identity;
    both are rejected here.  Returns (x, lcc).
    """
    from scipy.optimize import minimize
    from scipy.spatial.transform import Rotation

    def neg_obj(x):
        try:
            R_c = Rotation.from_rotvec(x[:3]).as_matrix()
        except ValueError:
            return 1.0
        if any(_rot_angle_between(R_c, S) < np.deg2rad(10.0)
               for S in sg_cart):
            return 1.0
        R_f = cart_to_frac(R_c, cell)
        t = x[3:6] - np.floor(x[3:6])
        cc, _ = _mean_local_cc(rho_d, core_d, R_f, t, win)
        return -cc

    res = minimize(neg_obj, np.asarray(x0, dtype=np.float64),
                   method='Nelder-Mead',
                   options=dict(maxiter=maxiter, xatol=1e-3, fatol=1e-3))
    cc = -float(res.fun)
    if cc <= 0.0:
        return None, -1.0
    return res.x, cc


def detect_ncs(rho, mask, cell, spacegroup, log=None, cc_min=0.15,
               core_frac=0.22, n_cand=10):
    """Detect a proper-rotation NCS operator from the current density.

    A plain protein mask covers roughly half the cell, and the overlap
    of such a mask with a *randomly* rotated copy of itself already
    scores ~0.5 -- there is no contrast in that score.  The search is
    therefore run on a tight mask of the molecular cores (top
    ``core_frac`` of the smoothed density), where the random baseline is
    ``core_frac`` and a true operator scores far above it.  The best
    rotation candidates are then refined (rotation + translation) by
    correlation against the downsampled map itself, and accepted only
    if the local correlation between the map and its averaged image --
    the quantity the averager actually relies on -- exceeds ``cc_min``.

    Returns (R_f, t_f, cc) or None.
    """
    from scipy.spatial.transform import Rotation

    rho = np.asarray(rho, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool)
    if log is None:
        log = (lambda *_: None)

    # --- core mask: top core_frac of the ~5 A smoothed density --------
    v = _voxel_size(cell, rho.shape)
    sm = ndimage.uniform_filter(rho, size=[max(3, int(round(5.0 / g)))
                                           for g in v], mode='wrap')
    core = (sm > np.percentile(sm, 100.0 * (1.0 - core_frac))) & mask
    if float(core.sum()) < 200:
        log('ncs: core mask too small to search')
        return None

    # --- downsampled grids --------------------------------------------
    f = max(1, int(round(min(rho.shape) / 24.0)))
    if f > 1:
        core_d = ndimage.zoom(core.astype(np.float64), 1.0 / f,
                               order=1) > 0.5
        rho_d = ndimage.zoom(rho, 1.0 / f, order=1)
    else:
        core_d, rho_d = core, rho
    cdf = core_d.astype(np.float64)
    if float(cdf.sum()) < 50:
        log('ncs: downsampled core mask too small')
        return None
    fr = max(1, f // 2)                    # refinement grid: 2x finer
    if fr > 1:
        rho_r = ndimage.zoom(rho, 1.0 / fr, order=1)
        core_r = ndimage.zoom(core.astype(np.float64), 1.0 / fr,
                               order=1) > 0.5
    else:
        rho_r, core_r = rho, core
    sg_rots = _sg_rotations(spacegroup)
    orth = np.asarray(cell.orth.mat.tolist(), dtype=np.float64)
    frac = np.asarray(cell.frac.mat.tolist(), dtype=np.float64)
    sg_cart = [orth @ r @ frac for r in sg_rots]

    def overlap_score(R_f, arr):
        """Best-translation normalised overlap of the core mask."""
        M, _ = _idx_matrix(R_f, np.zeros(3), arr.shape)
        rot = ndimage.affine_transform(
            arr, M, offset=np.zeros(3), order=1, mode='grid-wrap')
        c, _ = _fft_corr_peak(arr, rot)
        nrm = np.sqrt(float(np.sum(arr * arr))
                      * float(np.sum(rot * rot)) + 1e-12)
        return c / nrm

    # --- coarse scan (20 deg) on the core mask -------------------------
    step = np.deg2rad(20.0)
    cands = []
    for a in np.arange(0, 2 * np.pi, step):
        for b in np.arange(0, np.pi, step):
            for g in np.arange(0, 2 * np.pi, step):
                R_c = rot_zxz(a, b, g)
                R_f = cart_to_frac(R_c, cell)
                if any(np.allclose(R_f, r, atol=0.05) for r in sg_rots):
                    continue
                s = overlap_score(R_f, cdf)
                if s > core_frac * 1.15:   # just above the random baseline;
                    cands.append((s, R_c))  # multi-copy NCS only overlaps
                    # ~1/n_copies of the core, so the bar must stay low here
                    # and let the local-cc verification below do the work
    if not cands:
        log('ncs: no rotation candidate above the core-overlap baseline')
        return None
    cands.sort(key=lambda c: -c[0])
    log(f'ncs: {len(cands)} rotation candidates above baseline')

    # keep candidates whose rotations are well separated
    picked = []
    for s, Rc in cands:
        if any(_rot_angle_between(Rc, Rp) < np.deg2rad(20.0)
               for _, Rp in picked):
            continue
        picked.append((s, Rc))
        if len(picked) >= n_cand:
            break

    # windows for the local-correlation objective (in grid voxels)
    win_f = [max(2, min(15, int(round(8.0 / (g * f))))) for g in v]
    win_r = [max(2, min(15, int(round(8.0 / (g * fr))))) for g in v]
    win_full = [max(3, min(15, int(round(8.0 / g)))) for g in v]

    def polish(rho_d, core_d, R_c, win, maxiter):
        """6-parameter Nelder-Mead refinement of one candidate."""
        Rr = cart_to_frac(R_c, cell)
        M, _ = _idx_matrix(Rr, np.zeros(3), rho_d.shape)
        rot_r = ndimage.affine_transform(
            rho_d, M, offset=np.zeros(3), order=1, mode='grid-wrap')
        _, t_idx = _fft_corr_peak(rho_d, rot_r)
        t0 = np.asarray(t_idx, dtype=np.float64) / np.asarray(
            rho_d.shape, dtype=np.float64)
        x0 = np.concatenate([Rotation.from_matrix(R_c).as_rotvec(), t0])
        x, cc = _nm_refine(rho_d, core_d, cell, x0, sg_cart, win,
                           maxiter=maxiter)
        if x is None:
            return None
        R_c = Rotation.from_rotvec(x[:3]).as_matrix()
        if any(_rot_angle_between(R_c, S) < np.deg2rad(10.0)
               for S in sg_cart):
            return None
        return cc, cart_to_frac(R_c, cell), x[3:6] - np.floor(x[3:6])

    best = None                       # (cc, R_f, t_f)
    for s0, Rc in picked:
        # fine rotation scan (5 deg, rotvec space) around the coarse peak
        w0 = Rotation.from_matrix(Rc).as_rotvec()
        st = np.deg2rad(5.0)
        Rc_best, s_best = Rc, s0
        for ex in (-2 * st, -st, 0, st, 2 * st):
            for ey in (-2 * st, -st, 0, st, 2 * st):
                for ez in (-2 * st, -st, 0, st, 2 * st):
                    w2 = w0 + np.array([ex, ey, ez])
                    R2c = Rotation.from_rotvec(w2).as_matrix()
                    if any(_rot_angle_between(R2c, S) < np.deg2rad(10.0)
                           for S in sg_cart):
                        continue
                    s2 = overlap_score(cart_to_frac(R2c, cell), cdf)
                    if s2 > s_best:
                        Rc_best, s_best = R2c, s2
        # stage 1: coarse-grid local-correlation refinement
        out = polish(rho_d, core_d, Rc_best, win_f, maxiter=300)
        if out is None:
            continue
        cc1, R_f1, t_f1 = out
        log(f'ncs: candidate core-overlap {s0:.2f} -> local cc '
            f'{cc1:.3f} (coarse grid)')
        if cc1 < cc_min * 0.6:
            continue
        # stage 2: refinement-grid polish, then full-map verification
        orth = np.asarray(cell.orth.mat.tolist(), dtype=np.float64)
        frac_m = np.asarray(cell.frac.mat.tolist(), dtype=np.float64)
        R_c1 = orth @ R_f1 @ frac_m
        out = polish(rho_r, core_r, R_c1, win_r, maxiter=250)
        if out is not None:
            cc2, R_f2, t_f2 = out
        else:
            cc2, R_f2, t_f2 = cc1, R_f1, t_f1
        cc, reg = _mean_local_cc(rho, core, R_f2, t_f2, win_full,
                                 variants=True)
        log(f'ncs: candidate refined -> local cc {cc:.3f} (full map, '
            f'region {float(reg.mean()):.3f})')
        if cc >= cc_min and (best is None or cc > best[0]):
            best = (cc, R_f2, t_f2)
    if best is None:
        log(f'ncs: no operator with local cc >= {cc_min:.2f}; '
            'treating this crystal as non-NCS')
        return None
    cc, R_f, t_f = best
    log(f'ncs: detected operator local cc {cc:.3f}, translation '
        f'({t_f[0]:.3f}, {t_f[1]:.3f}, {t_f[2]:.3f}) frac')
    return R_f, t_f, cc
