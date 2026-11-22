#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""oasis_dm v2 -- the engine the CLI runs: cross-validated map + free-set stopping.

Why v2 exists
-------------
In v1 the Sim/sigmaA missing-structure variance SigmaQ was fitted on WORK
reflections, i.e. exactly those reflections used to compute the map.  That makes
(Fo - D*Fc) systematically too small, X too large, and the modified phases
over-trusted (FOM inflation).  v2 adds the two things that can be checked without
reference phases:

  * a second map per cycle built from the SAME mask and histogram LUT with the
    free reflections deleted, which gives a free-set amplitude CC / R (the dm
    COMBINE OMIT idea).  It is a *diagnostic*: the omit map still contains ~95% of
    the reflections, so its amplitude agreement is not an unbiased phase-error
    estimate, and `sigmaq_source='omit'` is consequently not the default.
  * free-set driven stopping: an absolute-tolerance plateau test with arg-max
    restore (dm NCYCLE AUTO / RESTORE idea).

Measured on the bundled cases: the free-set CC creeps up ~0.005 per cycle for as
long as you let it run, so a strict "no improvement" test never fires; an absolute
tolerance is required.  True phase error is flat from cycle 4-6 onwards.

Solvent masking, histogram matching, gamma correction, HL recombination and the
output columns are shared with v1 through ``core.DensityModifier``.  This program
never calls dm or cparrot.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from dataclasses import dataclass, asdict

import numpy as np

try:
    from . import __version__
    from .core import (Config, DensityModifier, hl_centroid, wrap180,
                       fom_to_x, x_to_fom, build_parser)
    from . import ncs
except ImportError:  # direct script execution
    from core import (Config, DensityModifier, hl_centroid, wrap180,
                      fom_to_x, x_to_fom, build_parser)
    import ncs


@dataclass
class ConfigV2(Config):
    sigmaq_source: str = 'blend'    # 'blend' (work-dominated, stable, default) |
                                    # 'kfold' (k-fold omit maps: every reflection
                                    # predicted from a map that never saw it --
                                    # the honest SigmaQ) |
                                    # 'omit' (unbiased attempt, needs scale
                                    # correction to be stable) | 'work'
    kfold: int = 5                  # number of k-fold omit maps per cycle
    extend: bool = True             # activate unphased amplitudes (phase
                                    # extension) progressively
    ext_steps: int = 5              # cycles over which extension activates
    aniso: bool = False             # per-cycle anisotropic Fo -> |Fmod| scaling
    ncs: str = ''                   # path to an NCS operator file (fractional)
    ncs_auto: bool = False          # detect the operator from the cycle-2 map,
                                    # gated on the free-set score
    auto_stop: bool = True
    max_cycles: int = 8
    min_cycles: int = 2
    patience: int = 2
    min_improve: float = 0.02       # free-set score gain that counts as progress
    free_metric: str = 'cc'         # 'cc' | 'rfree'
    min_free: float = 500.0         # below this many valid free reflections even
                                    # an explicit 'omit' request falls back to
                                    # blend (work-dominated SigmaQ)
    omit_fc: str = 'full'           # 'full' | 'omit': which |Fmod| predicts the
                                    # omitted reflections in the SigmaQ fit
    cmbwt: float = 0.4              # v1-calibrated damping; kept conservative
                                    # because the simplified omit SigmaQ estimate
                                    # is optimistic (free set is only ~5% of
                                    # reflections) and must not inflate FOM
    map_update: str = 'posterior'   # 'difference' iterates complex map coefficients
    map_output: str = 'standard'    # 'posterior' writes observed mFo and posterior phase
    n_shells: int = 10              # fewer shells -> more reflections each

    def as_dict(self):
        return asdict(self)


class DensityModifierV2(DensityModifier):
    """v1 engine + omit-based error estimation + free-set early stopping."""

    def _prepare(self):
        cfg = self.cfg
        if cfg.extend and np.any(self.ext):
            # the FFT grid must sample the finest *measured* amplitude, or the
            # extension predictions alias; with extension off, v2 keeps v1's
            # behaviour of sizing the grid to the phased data only
            self.dmin = min(self.dmin, float(np.nanmin(self.d[self.ext])))
        super()._prepare()
        self._active = self.valid     # phased + activated extension (run())
        self.free_valid = self.free & self.valid
        self.n_free_valid = int(self.free_valid.sum())
        if self.n_free_valid < 50:
            self.log(f'[warn] only {self.n_free_valid} valid free reflections; '
                     f'omit-based SigmaQ may be noisy')
        # Robustness guard: a handful of free reflections cannot pin down D and
        # SigmaQ per resolution shell.  Fall back to a work-dominated blend so we
        # do not inflate FOM (the symmetric failure mode to v1's self-reference
        # bias).  This keeps v2 honest on small-test-set SAD data.
        if self.cfg.sigmaq_source == 'omit' and self.n_free_valid < self.cfg.min_free:
            self.cfg.sigmaq_source = 'blend'
            self.cfg.cmbwt = 0.4
            self.log(f'[v2  ] free reflections {self.n_free_valid} < min_free '
                     f'{self.cfg.min_free}: SigmaQ source -> blend '
                     f'(work-dominated); cmbwt -> 0.4 to avoid FOM inflation')
        self.log(f'[v2  ] sigmaQ source={self.cfg.sigmaq_source}  '
                 f'auto_stop={self.cfg.auto_stop}  '
                 f'free reflections={self.n_free_valid}'
                 + (f'  kfold={self.cfg.kfold}' if self.cfg.sigmaq_source == 'kfold' else ''))

    # ------------------------------------------------------------------
    def _estimate_gamma(self, rho, modify):
        rng = np.random.default_rng(self.cfg.seed + 9000 + len(self.stats))
        eps = 0.03 * max(float(rho.std()), 1e-6)
        delta = rng.normal(0.0, eps, rho.shape)
        base = modify(rho)
        pert = modify(rho + delta)
        d = pert - base
        g = float(np.sum(d * delta) / max(float(np.sum(delta * delta)), 1e-12))
        return float(np.clip(g, 0.0, 0.9)), base

    # ------------------------------------------------------------------
    def _kfold_folds(self, k):
        """Deterministic, resolution-shell-balanced fold labels for work reflections.

        ``fold[i] in [0, k)`` for work reflections, ``-1`` otherwise.  Balanced
        per shell so every fold map misses the same fraction of every shell.
        """
        if getattr(self, '_fold', None) is not None and self._fold_k == k:
            return self._fold
        rng = np.random.default_rng(self.cfg.seed + 777)
        fold = np.full(self.nref, -1, dtype=np.int32)
        for s in range(self.cfg.n_shells):
            idx = np.nonzero(self.work & (self.shell == s))[0]
            perm = rng.permutation(len(idx))
            fold[idx[perm]] = np.arange(len(idx)) % k
        self._fold, self._fold_k = fold, k
        return fold

    def _map_with_deleted(self, F_map, phi_deg, modify, gamma):
        """Map coefficients -> (deleted set already zeroed) -> modified Fmod."""
        rho = self.rho_from_grid(self._grid_from_coeffs(F_map, phi_deg))
        rho_mod = modify(rho)
        if gamma is not None:
            rho_mod = (rho_mod - gamma * rho) / (1.0 - gamma)
        return self.fmod_from_rho(rho_mod)

    def _kfold_predictions(self, F_best, phi_deg, modify, gamma):
        """Predict every work reflection from a map that never contained it.

        The work reflections are split into ``k`` shell-balanced folds; fold
        ``j`` is predicted from the map built with fold ``j`` *and the free set*
        deleted.  The (Fo, Fpred) pairs therefore carry no self-reference bias,
        unlike a fit on the map's own reflections (v1) or on a ~5% free set
        predicted by a map that still contains the other ~95% of the data
        (``sigmaq_source='omit'``).  The per-shell D and SigmaQ fitted on them
        -- and the FOM that follows -- are honest.  Cost: ``k`` extra maps per
        cycle.

        Returns complex predictions; ``sigmaq_source='kfold'`` uses the
        magnitudes, ``'kcal'`` also the phases.
        """
        k = max(int(self.cfg.kfold), 2)
        fold = self._kfold_folds(k)
        Fpred = np.full(self.nref, np.nan, dtype=np.complex128)
        for j in range(k):
            sel = fold == j
            if not np.any(sel):
                continue
            F_map = F_best.copy()
            F_map[sel | self.free] = 0.0
            Fpred[sel] = self._map_with_deleted(F_map, phi_deg, modify, gamma)[sel]
        return Fpred

    # ------------------------------------------------------------------
    def _sim_weight_v2(self, Fc_full, Fc_omit, Fpred=None,
                       phi_best_rad=None, m_best=None):
        """Fit D and SigmaQ on OMITTED (free) reflections; returns X for all."""
        cfg = self.cfg
        X = np.zeros(self.nref, dtype=np.float64)
        Ds = np.zeros(cfg.n_shells)
        Qs = np.zeros(cfg.n_shells)
        Fo = np.nan_to_num(self.Fo)
        if cfg.sigmaq_source == 'work':
            fitmask = self.work
            fc_fit = Fc_full
        elif cfg.sigmaq_source in ('kfold', 'kcal') and Fpred is not None:
            # honest: every fitting reflection was predicted by a map built
            # without it, so (Fo - D*Fpred)^2 measures real map error
            ap = np.abs(Fpred)
            fitmask = self.work & np.isfinite(ap) & (ap > 0)
            fc_fit = ap
        elif cfg.sigmaq_source == 'blend':
            # work-dominated, free reflections only help the per-shell counts
            fitmask = self.work | self.free_valid
            fc_fit = Fc_full
        else:  # 'omit' (unbiased): fit on FREE reflections.  Fc may come from the
               # full map ('full') or from the map that never saw them ('omit').
               # The omit variant is the physically consistent one: a real unknown
               # reflection is absent from the map that predicts it.
            fitmask = self.free_valid
            fc_fit = Fc_omit if cfg.omit_fc == 'omit' else Fc_full
        for s in range(cfg.n_shells):
            idx = (fitmask & (self.shell == s) & np.isfinite(fc_fit)
                   & (fc_fit > 0) & (Fo > 0))
            n = int(idx.sum())
            if n < 8:
                Ds[s], Qs[s] = 1.0, np.nan
                continue
            fo = Fo[idx]
            fc = fc_fit[idx]
            D = float(np.sum(fo * fc) / np.sum(fc * fc))
            resid = float(np.mean((fo - D * fc) ** 2))
            Qs[s] = max(resid, 0.10 * float(np.mean(fo ** 2)))
            Ds[s] = D
        good = np.isfinite(Qs) & (Qs > 0)
        Dg = float(np.median(Ds[good])) if np.any(good) else 1.0
        Qg = (float(np.median(Qs[good])) if np.any(good)
              else float(np.mean(Fo[self.work] ** 2)))
        for s in range(cfg.n_shells):
            if not (np.isfinite(Qs[s]) and Qs[s] > 0):
                Ds[s], Qs[s] = Dg, Qg
        for s in range(cfg.n_shells):
            idx = self._active & (self.shell == s)
            if np.any(idx):
                X[idx] = np.clip(2.0 * Fo[idx] * (Ds[s] * Fc_full[idx]) / Qs[s],
                                 0.0, cfg.x_max)
        # --- kcal: rescale the per-shell concentration X to the map's ------
        # --- *measured* phase error, estimated from the disagreement -------
        # --- between out-of-map k-fold predictions and the current phases --
        # For independent errors E[cos(dphi_obs)] = cos_prior * cos_map, so
        # cos_map = cos_obs / cos_prior and X_shell = fom_to_x(cos_map).
        self._kcal_scale = np.ones(cfg.n_shells)
        if (cfg.sigmaq_source == 'kcal' and Fpred is not None
                and phi_best_rad is not None and m_best is not None):
            ap = np.abs(Fpred)
            dphi = np.angle(Fpred) - phi_best_rad
            w = np.where(self.valid, np.asarray(m_best, float), 0.0)
            for s in range(cfg.n_shells):
                idx = (self.work & (self.shell == s) & np.isfinite(ap)
                       & (ap > 0) & (w > 0.05))
                if idx.sum() < 8:
                    continue
                sw = float(np.sum(w[idx]))
                cos_obs = float(np.sum(w[idx] * np.cos(dphi[idx])) / sw)
                cos_prior = float(np.sum(w[idx] * m_best[idx]) / sw)
                r = cos_obs / max(cos_prior, 1e-3)
                x_target = fom_to_x(np.clip(r, 0.0, 0.9995))
                x_cur = float(np.mean(X[idx]))
                if np.isfinite(x_target) and x_cur > 1e-6:
                    self._kcal_scale[s] = x_target / x_cur
            scale = self._kcal_scale[np.clip(self.shell, 0, cfg.n_shells - 1)]
            X = np.where(self.valid, np.clip(X * scale, 0.0, cfg.x_max), 0.0)
        return X, Ds, Qs, Dg

    def _free_metric(self, Fc_omit):
        sel = self.free_valid & np.isfinite(Fc_omit) & (Fc_omit > 0)
        if sel.sum() < 20:
            return float('nan'), float('nan')
        fo = np.nan_to_num(self.Fo)[sel]
        fc = Fc_omit[sel]
        cc = float(np.corrcoef(fo, fc)[0, 1])
        k = float(np.sum(fo * fc) / np.sum(fc * fc))
        rfree = float(np.sum(np.abs(fo - k * fc)) / np.sum(fo))
        return cc, rfree

    def _fit_aniso(self):
        """Anisotropic amplitude scaling from direction statistics.

        Fits  ln Fo - <ln Fo>_shell ~ c + h^T B h  on work reflections: the
        per-shell mean is the isotropic reference, so B captures the
        direction-dependent falloff of the measured amplitudes.  (Fitting Fo
        against |Fmod| instead would be circular -- the map is built from Fo
        and inherits the anisotropy, leaving no residual to detect.)

        Returns the 7 coefficients [c, B11, B22, B33, B12, B13, B23] or None
        when there is too little data.  h are integer Miller indices: the
        anisotropy is tied to the crystal axes, so fractional indices are the
        natural coordinate system.
        """
        cfg = self.cfg
        sel = self.work & np.isfinite(self._Fo_raw) & (self._Fo_raw > 0)
        if int(sel.sum()) < 200:
            return None
        h = self.hkl[sel].astype(np.float64)
        lfo = np.log(self._Fo_raw[sel])
        sh = self.shell[sel]
        mu = np.zeros(cfg.n_shells)
        for s in range(cfg.n_shells):
            m = sh == s
            if np.any(m):
                mu[s] = lfo[m].mean()
        y = lfo - mu[sh]
        A = np.column_stack([np.ones(len(h)),
                             h[:, 0] ** 2, h[:, 1] ** 2, h[:, 2] ** 2,
                             2 * h[:, 0] * h[:, 1], 2 * h[:, 0] * h[:, 2],
                             2 * h[:, 1] * h[:, 2]])
        # small ridge on the quadratic terms keeps B tame on thin shells
        ridge = np.zeros((7, 7))
        ridge[1:, 1:] = np.eye(6) * 1e-6 * len(h)
        coef = np.linalg.solve(A.T @ A + ridge, A.T @ y)
        return coef

    def _aniso_factor(self, coef):
        """exp(-(c + h^T B h)) per reflection, capped to [0.5, 2]."""
        h = self.hkl.astype(np.float64)
        quad = (coef[1] * h[:, 0] ** 2 + coef[2] * h[:, 1] ** 2
                + coef[3] * h[:, 2] ** 2 + coef[4] * 2 * h[:, 0] * h[:, 1]
                + coef[5] * 2 * h[:, 0] * h[:, 2] + coef[6] * 2 * h[:, 1] * h[:, 2])
        lnfac = -(coef[0] + quad)
        return np.exp(np.clip(lnfac, -np.log(2.0), np.log(2.0)))

    # ------------------------------------------------------------------
    def run(self):
        cfg = self.cfg
        HL_exp = self.HL_exp.copy()
        use_difference_map = cfg.map_update == 'difference'
        if use_difference_map and self.hl_source.startswith('PHIB/FOM'):
            centric = self.valid & self.centric
            field = np.arctanh(np.clip(self.fom_in[centric], 0, 0.99995))
            phase = self._project_centric(np.deg2rad(self.phi_exp_deg))[centric]
            HL_exp[centric, 0] = field * np.cos(phase)
            HL_exp[centric, 1] = field * np.sin(phase)
        HL_cur = HL_exp.copy()
        moments = self._phase_moments if use_difference_map else hl_centroid
        Fo = np.nan_to_num(self.Fo)
        stats = []
        self.stats = stats          # shared list, used for per-cycle RNG seeding
        self.early_stop_available = None
        rho_prev = None
        t0 = time.time()
        ncycles = cfg.max_cycles if cfg.auto_stop else cfg.cycles
        has_free_observations = bool(np.any(self.free_valid))
        best = None
        prev_score = None
        plateau = 0
        free_scored = False
        # --- optional NCS averaging --------------------------------------
        self._ncs_av = None
        self._ncs_tried = False
        if cfg.ncs:
            ops = ncs.read_operators(cfg.ncs, self.cell)
            self._ncs_av = ncs.NcsAverager(ops, spacegroup=self.sg, cell=self.cell)
            self.log(f'[v2  ] NCS averaging: {len(ops)} operator(s) from '
                     f'{cfg.ncs}')
        self._Fo_raw = self.Fo.copy()       # anisotropic correction is always
                                            # applied to the RAW amplitudes
        # --- optional anisotropic scaling of Fo (once; data-only estimate) --
        self._aniso = None
        if cfg.aniso:
            coef = self._fit_aniso()
            if coef is not None:
                fac = self._aniso_factor(coef)
                meas = np.isfinite(self._Fo_raw) & (self._Fo_raw > 0)
                self.Fo = np.where(meas, self._Fo_raw * fac, self.Fo)
                self._aniso = dict(B=[float(x) for x in coef],
                                   fac_min=float(fac[meas].min()),
                                   fac_max=float(fac[meas].max()))
                self.log(f'[v2  ] anisotropic correction: B_diag='
                         f'{coef[1]:+.3g} {coef[2]:+.3g} {coef[3]:+.3g}; '
                         f'Fo rescaled by {self._aniso["fac_min"]:.2f}-'
                         f'{self._aniso["fac_max"]:.2f} (capped at 0.5-2.0)')
            else:
                self.log('[v2  ] anisotropic correction: too few work '
                         'reflections; skipped')

        # --- phase extension: unphased amplitudes activate progressively ----
        ext_work = self.ext & ~self.free
        n_ext = int(ext_work.sum())
        ext_act = np.zeros(self.nref, dtype=bool)
        ext_order = None
        n_on_prev = 0
        previous_active = self.valid
        if n_ext:
            if cfg.extend:
                order = np.argsort(-self.d[ext_work])     # coarse shells first
                ext_order = np.nonzero(ext_work)[0][order]
                self.log(f'[v2  ] phase extension: {n_ext} reflections have Fo '
                         f'but no starting phase (d '
                         f'{np.nanmin(self.d[ext_work]):.2f}-'
                         f'{np.nanmax(self.d[ext_work]):.2f} A); activating the '
                         f'coarsest {100 // max(cfg.ext_steps, 1)}% per cycle '
                         f'over {cfg.ext_steps} cycles')
            else:
                self.log(f'[v2  ] {n_ext} unphased amplitudes present; phase '
                         f'extension disabled (--no-extend)')

        self.log('')
        self.log(f'{"cyc":>3} {"meanFOM":>8} {"freeCC":>7} {"freeR":>7} '
                 f'{"SigmaQ":>9} {"gamma":>7} {"dPhi":>7} {"mapCC":>7} {"best":>5}')
        for cyc in range(1, ncycles + 1):
            # --- phase-extension activation (coarse shells first) ----------
            if ext_order is not None:
                n_on = min(int(round(n_ext * min(1.0, cyc / max(cfg.ext_steps, 1)))),
                           n_ext)
                if n_on > n_on_prev:
                    ext_act[:] = False
                    ext_act[ext_order[:n_on]] = True
                    self.log(f'[v2  ] cycle {cyc}: activated {n_on - n_on_prev} '
                             f'more extension reflections ({n_on}/{n_ext} total, '
                             f'down to {self.d[ext_order[n_on - 1]]:.2f} A)')
                    n_on_prev = n_on
            active = self.valid | ext_act
            self._active = active
            Fo = np.nan_to_num(self.Fo)     # refreshed: aniso may rescale Fo

            phi_best_rad, m_best = moments(HL_cur)
            phi_best_rad = self._project_centric(phi_best_rad)
            phi_deg = np.rad2deg(phi_best_rad)
            F_best = np.where(active, m_best * Fo, 0.0)
            if use_difference_map and cyc > 1:
                best_map, _ = self._difference_coefficients(
                    F_best, phi_best_rad, Fmod, Ds[self.shell])
                best_map = np.where(previous_active, best_map, 0j)
                F_best, phi_deg = np.abs(best_map), np.rad2deg(np.angle(best_map))
            previous_active = active
            F_omit = np.where(self.free, 0.0, F_best) if has_free_observations else None

            # --- full map -------------------------------------------------
            rho = self.rho_from_grid(self._grid_from_coeffs(F_best, phi_deg))
            radius = max(cfg.radius_min, cfg.radius * cfg.radius_decay ** (cyc - 1))
            mean_fom_now = float(np.mean(m_best[self.valid]))
            modify, mask_solv, mask_prot, d_eff = self._build_modifier(
                rho, radius, mean_fom_now)
            # --- optional NCS averaging (wraps the modification) ---------
            # detection waits for cycle 2: one solvent-flattening round
            # first lifts the map SNR enough for the local-correlation
            # search to see the NCS mates
            if (cfg.ncs_auto and self._ncs_av is None
                    and not self._ncs_tried and cyc >= 2):
                self._ncs_tried = True
                det = ncs.detect_ncs(rho, mask_prot > 0.5, self.cell, self.sg,
                                    log=self.log)
                if det is not None:
                    R_f, t_f, score = det
                    av = ncs.NcsAverager(
                        [(R_f, t_f)], spacegroup=self.sg, cell=self.cell)
                    self.log(f'[v2  ] NCS auto-detection: operator found '
                             f'(score {score:.3f})')
                    # --- free-set gate -------------------------------------
                    # Averaging with a mediocre operator, or on a map that
                    # plain DM already solves, can hurt.  Trial both
                    # variants on THIS cycle's map and keep averaging
                    # only if it does not degrade the free-set cc.
                    def _trial_cc(mdf):
                        if F_omit is None:
                            return float('nan')
                        g = None
                        if cfg.gamma_on:
                            g, rmod = self._estimate_gamma(rho, mdf)
                        else:
                            rmod = mdf(rho)
                        if g is not None:
                            rmod = (rmod - g * rho) / (1.0 - g)
                        fmo = self._map_with_deleted(F_omit, phi_deg, mdf, g)
                        return self._free_metric(np.abs(fmo))[0]

                    cc_plain = _trial_cc(modify)
                    cc_ncs = _trial_cc(av.wrap(modify, mask_prot > 0.5))
                    # the amplitude-only free metric cannot see phase
                    # damage from a mediocre operator, so a bare >= is
                    # not enough: on 7tj1 a single half-detected
                    # operator scored +0.014 free cc and cost 0.07 cos
                    # of truth.  demand a clear margin.
                    margin = 0.03
                    if not (np.isfinite(cc_plain) and np.isfinite(cc_ncs)):
                        self._ncs_av = av   # no usable free set: trust NCS
                        self.log('[v2  ] NCS gate: no usable free set; '
                                 'keeping averaging')
                    elif cc_ncs > cc_plain + margin:
                        self._ncs_av = av
                        self.log(f'[v2  ] NCS gate: free cc {cc_ncs:.3f} '
                                 f'(with NCS) > {cc_plain:.3f} (plain) + '
                                 f'{margin}; keeping averaging')
                    else:
                        self.log(f'[v2  ] NCS gate: free cc {cc_ncs:.3f} '
                                 f'(with NCS) <= {cc_plain:.3f} (plain) + '
                                 f'{margin}; falling back to plain density '
                                 'modification')
                else:
                    self.log('[v2  ] NCS auto-detection: no operator found, '
                             'continuing without averaging')
            if self._ncs_av is not None:
                modify = self._ncs_av.wrap(modify, mask_prot > 0.5)
            gamma = None
            if cfg.gamma_on:
                gamma, rho_mod = self._estimate_gamma(rho, modify)
            else:
                rho_mod = modify(rho)
            if gamma is not None:
                rho_mod = (rho_mod - gamma * rho) / (1.0 - gamma)
            Fmod = self.fmod_from_rho(rho_mod)
            Fc_full = np.abs(Fmod)

            # --- omit map (SAME mask, SAME LUT) ---------------------------
            if F_omit is not None:
                Fc_omit = np.abs(self._map_with_deleted(F_omit, phi_deg, modify, gamma))
                free_cc, free_r = self._free_metric(Fc_omit)
            else:
                Fc_omit = Fc_full
                free_cc = free_r = float('nan')

            # --- k-fold omit maps: honest out-of-map predictions ----------
            Fpred = None
            if cfg.sigmaq_source in ('kfold', 'kcal'):
                Fpred = self._kfold_predictions(F_best, phi_deg, modify, gamma)

            # --- unbiased error estimate + phase recombination ------------
            X, Ds, Qs, Dg = self._sim_weight_v2(Fc_full, Fc_omit, Fpred,
                                                phi_best_rad, m_best)
            wgt = cfg.cmbwt
            phi_mod = np.angle(Fmod)
            HL_mod = np.column_stack([wgt * X * np.cos(phi_mod),
                                      wgt * X * np.sin(phi_mod),
                                      np.zeros(self.nref), np.zeros(self.nref)])
            HL_mod = np.where(active[:, None], HL_mod, 0.0)
            if cfg.anchor == 'exp':
                HL_new = HL_exp + HL_mod      # anchored to the ORIGINAL experiment
            else:
                HL_new = HL_cur + HL_mod      # cumulative (faster, more model bias)
            phi_new_rad, m_new = moments(HL_new)
            phi_new_rad = self._project_centric(phi_new_rad)
            HL_cur = HL_new

            # --- diagnostics ----------------------------------------------
            sel = self.valid
            dphi = float(np.mean(np.abs(wrap180(np.rad2deg(phi_new_rad[sel])
                                                - np.rad2deg(phi_best_rad[sel])))))
            map_cc = (float(np.corrcoef(rho_mod.ravel(), rho_prev.ravel())[0, 1])
                      if rho_prev is not None else float('nan'))
            rho_prev = rho_mod
            qmed_vals = Qs[np.isfinite(Qs) & (Qs > 0)]
            qmed = float(np.median(qmed_vals)) if qmed_vals.size else float('nan')
            score = free_cc if cfg.free_metric == 'cc' else -free_r
            snapshot = dict(HL=HL_new.copy(), phi=phi_new_rad.copy(), m=m_new.copy(),
                            Fmod=Fmod.copy(), Dg=Dg, Ds=np.array(Ds, copy=True),
                            Qs=np.array(Qs, copy=True),
                            X=X.copy(), ext_act=ext_act.copy(),
                            score=score, cycle=cyc)
            if not np.isfinite(score):
                # No usable free set: the early-stop criterion does not exist.
                # Stopping here would freeze the very first cycle, which is the
                # least converged iterate, so keep the latest one and run to the
                # requested cycle count instead.
                self.early_stop_available = False
                if not free_scored:
                    self.log('[v2  ] no valid free reflections: the early-stop '
                             'criterion is unavailable; running the full cycle '
                             'count and keeping the final iterate')
                    free_scored = True
                best = snapshot
                is_best = True
            elif best is None or score > best['score']:
                self.early_stop_available = True
                best = snapshot
                is_best = True
            else:
                is_best = False
            # Plateau test on the free-set score itself.  An absolute tolerance
            # is essential: the free-set CC of a density-modified map creeps up
            # by ~0.005 per cycle forever, so a strict "no improvement" test
            # never fires and the run always burns max_cycles.
            if np.isfinite(score):
                if prev_score is not None and score - prev_score < cfg.min_improve:
                    plateau += 1
                elif prev_score is not None:
                    plateau = 0
                prev_score = score

            stats.append(dict(cycle=cyc, mean_fom=float(np.mean(m_new[sel])),
                              mean_fom_ext=(float(np.mean(m_new[ext_act]))
                                            if np.any(ext_act) else None),
                              n_ext_active=int(ext_act.sum()),
                              aniso=self._aniso,
                              ncs_avg_frac=(self._ncs_av.last_frac
                                            if self._ncs_av is not None
                                            else None),
                              free_cc=free_cc, free_r_omit=free_r,
                              sigmaq_median=qmed, gamma=gamma, d_eff=d_eff,
                              mask_radius=radius,
                              solvent_frac=float(mask_solv.mean()),
                              map_cc_prev=map_cc, mean_dphi_deg=dphi,
                              mean_X=float(np.mean(X[sel])),
                              is_best=is_best, best_cycle=best['cycle'],
                              D_shells=[float(x) for x in Ds],
                              SigmaQ_shells=[float(x) for x in Qs],
                              shell_d=[float(x) for x in self.shell_d]))
            self.log(f'{cyc:3d} {np.mean(m_new[sel]):8.4f} {free_cc:7.4f} '
                     f'{free_r:7.4f} {qmed:9.2f} '
                     f'{(gamma if gamma is not None else float("nan")):7.3f} '
                     f'{dphi:7.2f} {map_cc:7.4f} {"*" if is_best else " ":>5}')

            ext_done = (ext_order is None) or (int(ext_act.sum()) >= n_ext)
            if (cfg.auto_stop and cyc >= cfg.min_cycles and plateau >= cfg.patience
                    and ext_done):
                self.log(f'[v2  ] early stop at cycle {cyc}: free-set '
                         f'{cfg.free_metric} improved by less than '
                         f'{cfg.min_improve} for {plateau} cycles; '
                         f'restoring cycle {best["cycle"]}')
                break

        # restore best cycle -------------------------------------------------
        self.best_phi_rad = best['phi']
        self.best_m = best['m']
        self.HL_final = best['HL']
        self.Fmod = best['Fmod']
        self.Dg = best['Dg']
        self.D_shell = best['Ds']
        self.stats = stats
        self.best_cycle = best['cycle']
        self.ext_activated = best['ext_act']
        # Extension reflections that never activated (the free set, or cycles
        # ran out): phase and FOM from the best map's own prediction.
        miss = self.ext & ~best['ext_act'] if cfg.extend else np.zeros(0, bool)
        if np.any(miss):
            phi_map = self._project_centric(np.angle(best['Fmod']))
            # X was only ever computed for active reflections; rebuild it for
            # the never-activated ones from the best cycle's shell D / SigmaQ
            Ds = np.asarray(best['Ds'])
            Qs = np.asarray(best['Qs'])
            sh = np.clip(self.shell[miss], 0, cfg.n_shells - 1)
            Fc_m = np.abs(best['Fmod'][miss])
            Fo_m = np.nan_to_num(self.Fo[miss])
            x = np.clip(cfg.cmbwt * 2.0 * Fo_m * Ds[sh] * Fc_m / Qs[sh],
                        0.0, cfg.x_max)
            self.best_phi_rad[miss] = phi_map[miss]
            self.best_m[miss] = np.clip(x_to_fom(x), 0.0, 1.0)
            HLx = np.zeros((int(miss.sum()), 4))
            HLx[:, 0] = x * np.cos(phi_map[miss])
            HLx[:, 1] = x * np.sin(phi_map[miss])
            self.HL_final[miss] = HLx
        self.elapsed = time.time() - t0
        self._write_output()
        self._write_stats()
        return self.stats

    def _write_output(self):
        # Extension reflections carry a measured Fo and now a phase, so they
        # must be written like any other reflection (FDM = Fo, filled
        # PHIDM/FOMDM) instead of falling through v1's "no data" branch.
        # With extension disabled, v1's behaviour is kept untouched.
        if self.cfg.extend and np.any(self.ext):
            orig = self.valid
            self.valid = self.valid | self.ext
            try:
                super()._write_output()
            finally:
                self.valid = orig
        else:
            super()._write_output()

    def _write_stats(self):
        # add v2-specific fields to the parent JSON payload
        super()._write_stats()
        if not self.cfg.stats or not os.path.isfile(self.cfg.stats):
            return
        import json
        with open(self.cfg.stats) as fh:
            payload = json.load(fh)
        payload['final'].update(dict(best_cycle=self.best_cycle,
                                     n_cycles_run=len(self.stats),
                                     sigmaq_source=self.cfg.sigmaq_source,
                                     map_update=self.cfg.map_update,
                                     map_output=self.cfg.map_output,
                                     auto_stop=self.cfg.auto_stop,
                                     early_stop_available=self.early_stop_available,
                                     n_ext=int(np.sum(self.ext)),
                                     n_ext_activated=int(np.sum(self.ext_activated)),
                                     aniso=self.cfg.aniso,
                                     ncs=(self.cfg.ncs or None),
                                     ncs_auto=self.cfg.ncs_auto,
                                     best_free_cc=max(
                                         (s['free_cc'] for s in self.stats
                                          if np.isfinite(s['free_cc'])),
                                         default=float('nan'))))
        with open(self.cfg.stats, 'w') as fh:
            json.dump(payload, fh, indent=2, ensure_ascii=False)


def build_parser_v2():
    p = build_parser()
    p.description = ('oasis_dm v2: density modification with omit-based SigmaQ '
                     'estimation and free-set early stopping')
    p.add_argument('--version', action='version',
                   version='oasis-dm ' + __version__)
    p.add_argument('--sigmaq-source',
                   choices=['kfold', 'kcal', 'omit', 'work', 'blend'],
                   default='blend',
                   help='which reflections fit SigmaQ (default blend = '
                        'work-dominated, stable; kfold = k-fold omit maps, '
                        'every reflection predicted from a map that never saw '
                        'it, giving an honest SigmaQ and FOM; kcal = kfold '
                        'plus a phase-disagreement calibration of the X '
                        'concentration; omit is an experimental unbiased '
                        'attempt)')
    p.add_argument('--map-update', choices=['posterior', 'difference'],
                   default='posterior',
                   help='Recycle posterior mFo or complex 2mFo-DFc coefficients')
    p.add_argument('--map-output', choices=['standard', 'posterior'],
                   default='standard',
                   help='Write standard map coefficients or posterior mFo/phase in DM_FWT/DM_PHWT')
    p.add_argument('--kfold', type=int, default=5,
                   help='number of k-fold omit maps per cycle (kfold mode)')
    p.add_argument('--no-extend', dest='extend', action='store_false',
                   default=True,
                   help='disable phase extension of reflections that have Fo '
                        'but no starting phase (NaN PHIB/FOM); by default they '
                        'are activated progressively, coarse shells first')
    p.add_argument('--ext-steps', type=int, default=5,
                   help='cycles over which extension reflections activate')
    p.add_argument('--aniso', action='store_true', default=False,
                   help='per-cycle anisotropic scaling of Fo towards |Fmod| '
                        '(6-parameter fit; the corrected amplitudes are used '
                        'internally and written to FDM)')
    p.add_argument('--ncs', default='',
                   help='NCS operator file: per operator 3 lines of '
                        '"r11 r12 r13 t1" (fractional rotation row + '
                        'translation); averaging is applied inside the '
                        'protein-mask intersection')
    p.add_argument('--ncs-auto', action='store_true', default=False,
                   help='detect an NCS operator automatically from the '
                        'cycle-2 map (SO(3) scan + FFT translation search) '
                        'and keep averaging only if the free-set score '
                        'beats plain density modification by a clear '
                        'margin (amplitude-only scores cannot see phase '
                        'damage, so ties fall back to plain)')
    p.add_argument('--auto-stop', action='store_true', default=True)
    p.add_argument('--no-auto-stop', dest='auto_stop', action='store_false',
                   help='disable free-set early stopping (fixed --cycles)')
    p.add_argument('--max-cycles', type=int, default=8)
    p.add_argument('--min-cycles', type=int, default=2)
    p.add_argument('--patience', type=int, default=2)
    p.add_argument('--free-metric', choices=['cc', 'rfree'], default='cc')
    p.add_argument('--hist-mode', choices=['cell', 'protein'], default='cell',
                   help='histogram target: whole-cell hist.lib distribution or the '
                        'protein-only part of it (default cell)')
    p.add_argument('--per-shell-d', action='store_true',
                   help='use the per-resolution-shell sigma_A in FWT/FDM output')
    p.add_argument('--omit-fc', choices=['full', 'omit'], default='full',
                   help='which |Fmod| predicts omitted reflections in the SigmaQ fit')
    p.add_argument('--min-improve', type=float, default=0.02,
                   help='free-set score gain that counts as progress (plateau test)')
    p.add_argument('--min-free', type=float, default=500.0,
                   help='below this many valid free reflections, auto-fallback '
                        'SigmaQ to blend (work-dominated) + cmbwt 0.4')
    p.set_defaults(n_shells=10)
    return p


def main(argv=None):
    args = build_parser_v2().parse_args(argv)
    radius = args.radius if args.radius is not None else (8.0 if args.mask == 'wang' else 4.0)
    radius_min = (args.radius_min if args.radius_min is not None
                  else (4.0 if args.mask == 'wang' else 3.0))
    cfg = ConfigV2(mtzin=args.mtzin, mtzout=args.mtzout,
                   fo=args.fo, sigf=args.sigf, phib=args.phib, fom=args.fom,
                   hl=args.hl, free=args.free, free_marker=args.free_marker,
                   truth=args.truth,
                   solc=args.solc, cycles=args.cycles, mask=args.mask,
                   radius=radius, radius_min=radius_min,
                   radius_decay=args.radius_decay, mask_smooth=args.mask_smooth,
                   hist_lib=args.hist_lib, hist_on=not args.no_hist,
                   solv_on=not args.no_solv, squash=not args.no_squash,
                   gamma_on=not args.no_gamma, omit_free=args.omit_free,
                   anchor=args.anchor, n_shells=args.n_shells, x_max=args.x_max,
                   cmbwt=(args.cmbwt if args.cmbwt is not None else 0.4),
                   restore=args.restore,
                   sample_rate=args.sample_rate, seed=args.seed,
                   stats=args.stats, tag=args.tag, quiet=args.quiet,
                   sigmaq_source=args.sigmaq_source, auto_stop=args.auto_stop,
                   max_cycles=args.max_cycles, min_cycles=args.min_cycles,
                   patience=args.patience, free_metric=args.free_metric,
                   min_free=args.min_free, min_improve=args.min_improve,
                   omit_fc=args.omit_fc, hist_mode=args.hist_mode,
                   per_shell_d=args.per_shell_d, kfold=args.kfold,
                   extend=args.extend, ext_steps=args.ext_steps,
                   aniso=args.aniso, ncs=args.ncs, ncs_auto=args.ncs_auto,
                   map_update=args.map_update, map_output=args.map_output)
    DensityModifierV2(cfg).run()
    return 0


if __name__ == '__main__':
    sys.exit(main())
