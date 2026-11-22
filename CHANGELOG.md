# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[semantic versioning](https://semver.org/).

## [2.4.0] — 2026-11-23

A real desktop window — the browser GUI from 2.3.0 is withdrawn.

### Changed

- **`oasis-dm-gui` is now a native desktop application** (PySide6 / Qt 6)
  instead of a browser page. Double-clicking the window is the normal way in;
  it does not need a browser, a server, or a terminal. The 2.3.0 browser front
  end was short-lived: it worked, but "a web page on localhost" is not what
  people mean by a GUI, and it made the run history awkward to inspect. The
  window keeps the same contract — one control per real CLI flag, only
  non-default options emitted, the exact command line previewed live, the
  engine log streamed, free-set cc plotted per cycle, Stop kills the engine
  cleanly — and adds a native file picker and column auto-detection with a
  header summary (`name · spacegroup · N reflections · d_min`).

### Added

- **`tools/make_macos_app.py`** — builds `dist/OasisDM.app`, a self-contained
  double-clickable bundle that embeds its own virtual environment with PySide6
  and the package, so it runs on a Mac without any Python setup. Verified by
  launching the built bundle. Not signed; fine for local or lab use, and the
  script says so.
- Optional extra `pip install ".[gui]"` (PySide6-Essentials ≥ 6.5,
  LGPL — the "Essentials" wheel, not the full PySide6, which drags in
  WebEngine and triples the bundle size). The CLI and the
  engine never import it, and the test suite skips the window tests when it is
  missing, so CI stays toolkit-free.

### Removed

- `oasis_dm/gui_page.py` and the local HTTP server (2.3.0's browser interface).


## [2.3.0] — 2026-11-23

A graphical interface, plus a version-string fix.

### Added

- **`oasis-dm-gui` (or `python -m oasis_dm.gui`)** — a browser-based front end.
  It starts a local HTTP server from the standard library only (bound to
  127.0.0.1, per-session token, no external assets or JavaScript libraries)
  and opens the default browser. The page offers a directory browser for
  picking the input/output MTZ, auto-detects the Fo / σ(Fo) / PHIB / FOM /
  free columns from the file header, exposes the engine options (solvent
  content, max cycles, free-set early stop, mask, histogram target, free
  metric, phase extension, anisotropic correction, NCS none/file/auto,
  SigmaQ source, histogram library), previews the exact `oasis-dm` command
  line before launching, streams the engine log, plots free-set cc per cycle,
  and stops the run cleanly. The engine always runs as a
  `python -m oasis_dm` subprocess, so the GUI adds no dependency and no GUI
  toolkit on the Python side — deliberately not Tkinter, whose system-Python
  builds are frequently unusable (on macOS 15 with the Xcode toolchain
  Python, `Tk()` aborts the process).
  - Tests: argv construction and column auto-detection logic, engine log
    parsing, directory listing, and full HTTP round trips including a real
    engine run and a mid-run Stop through the same JSON API the page uses.

### Fixed

- `oasis_dm.__version__` was still `2.1.0` while the package metadata said
  2.2.0, so `oasis-dm --version` reported the wrong release. Both now say
  2.3.0.

## [2.2.0] — 2026-11-22

Three engine capabilities: honest k-fold FOM calibration, phase extension with
anisotropic correction, and NCS averaging. Everything below was measured
against reference phases on the bundled cases.

### Added

- **NCS averaging (`--ncs FILE`).** Operator files list proper rotations in
  fractional coordinates (`x' = R_f x + t_f`, three `r r r t` rows per
  operator); improper rotations are rejected. The operator set is closed under
  space-group conjugation and inverse walks, so open A→B-only files work and
  the averaged map stays space-group symmetric. Averaging runs after solvent
  flattening and histogram matching inside every cycle; the perturbation-γ
  estimate and the k-fold Σ_Q predictions both see the averaged map, so the
  honesty properties below carry over.
  - Operators fitted on unwrapped coordinates obey `O(f + L) = O(f) + R·L`;
    all 27 lattice placements `t + R·L` are evaluated per operator power (the
    grid cannot be pre-rolled — `R·L` is not an integer offset in general).
  - Mask regions are partitioned by Voronoi growth from confident seeds
    (double-smoothed local correlation > 0.30), following the molecular copies
    even when the protein mask is one connected blob.
  - The blend weight is the double-smoothed local correlation of the assigned
    variant (≤ 0.10 keeps the un-averaged density, ≥ 0.25 replaces it fully):
    a failed partner search degrades that region only instead of smearing the
    map.
- **`--ncs-auto`.** From cycle 2 the engine can detect a proper-rotation
  operator from the current map (core-mask rotation scan + Nelder–Mead
  rotation/translation refinement) and enables it only when the free-set
  correlation beats plain density modification by a clear margin (+0.03;
  amplitude-only free scores cannot see phase damage, so ties fall back to
  plain). A bare >= comparison was tried and rejected: on `7tj1` it kept a
  single half-detected operator (local cc 0.32, one of the three true ones)
  that looked +0.014 better on free amplitudes yet cost 0.07 cos of truth
  (0.465 vs 0.538 plain). With the margin the gate declines to average on
  every bundled case — `7e1d` (spurious operator), `3i0t_hard`, `7tj1` —
  making it a conservative safety net rather than a substitute for the
  expert operator file.
- **k-fold omit maps (`--sigmaq-source kfold`).** Every work reflection is
  predicted from a map built without it (shell-balanced folds, free set always
  deleted). Fitting D and Σ_Q on these out-of-map predictions removes the
  self-reference bias that inflated FOM in v1. Measured on the bundled cases
  (cmbwt 0.4): FOM optimism 0.086/0.253/0.128 vs 0.106/0.296/0.145 for the
  blend estimate, at the same phase accuracy. `--kcal` additionally rescales
  the per-shell von Mises concentration from the k-fold phase disagreement;
  it helps `6zpp` but not the other cases, so it stays opt-in/experimental.
- **Phase extension (on by default, `--no-extend` to disable).** Reflections
  with `Fo` but no starting phase (NaN `PHIB`/`FOM`) are activated
  progressively, coarse shells first over `--ext-steps` cycles (default 5);
  never-activated extension reflections in the free set get phase/FOM from
  the best map prediction. Early stopping waits for full activation.
- **`--aniso` — anisotropic amplitude correction.** One-shot 6-parameter
  scaling of `Fo` from direction statistics (`ln Fo − ⟨ln Fo⟩_shell ∼ c +
  h'Bh`), capped at 0.5–2.0×. (Fitting `Fo` against `|Fmod|` was tried and
  rejected as circular — the map inherits its own anisotropy.)
- **NCS test cases** `3i0t`/`3i0t_hard` (2-fold) and `7tj1` (4-copy, 3
  operators) with operator files, and `tools/make_test_case.py --ncs-out` to
  regenerate them or build new ones from any PDB entry; `tests/test_ncs.py`.

### Measured (NCS)

Mean cos phase error vs reference (`PHIC`), work reflections:

| case | plain | `--ncs` |
| --- | ---: | ---: |
| 7tj1 (4-copy, operators 0.16–0.46 Å rmsd, wMPE 58°) | 0.538 | 0.552 |
| 3i0t (2-fold, operator 1.02 Å rmsd, wMPE 53°) | 0.728 | 0.679 |
| 3i0t_hard (2-fold, operator 1.02 Å rmsd, wMPE 70°) | 0.358 | 0.330 |

NCS averaging pays when the operators are tight and the copies are many (7tj1:
+0.014 cos); a loose operator costs accuracy and the poorer the starting phases
the deeper the loss (3i0t/3i0t_hard). That asymmetry is why `--ncs` is
documented as an expert option and `--ncs-auto` gates on the free set with a
clear margin.

### Fixed

- First run of the test suite on Windows/CI: `test_missing_hist_lib` left the
  process inside its temp directory (Windows refuses to delete the cwd), and
  `test_extension_beats_no_extend` relied on auto-stop cycle choice, which
  varies with the platform's FFT numerics. The extension effect is now
  measured at a pinned cycle count (+0.012 cos, identical on macOS, Linux and
  Windows) with a threshold to match.

## [2.1.0] — 2026-11-21

Focused on correctness, measured rather than assumed. Everything below was
checked against reference phases on the three bundled cases.

### Fixed

- **Bundled `6zpp.mtz` had an inconsistent truth column.** `PHIC` was related to
  `PHIB` by a 0.5 origin shift along **b**, which made the case look as if the
  engine destroyed the phases (weighted mean phase error 92°, map CC −0.02).
  Corrected (`PHIC -> PHIC + 180·k`); the starting phases now score 59.5° / 0.298.
- **`--anchor` was ignored by the v2 engine.** The option was accepted and then
  hard-coded to the `exp` branch, so `--anchor cumulative` silently did nothing
  (and would have inflated the FOM to ~0.88 if it had).
- **Duplicate `DM_*` columns when re-running on an earlier output.** The output
  block now rebuilds the file and drops stale `DM_*` labels by name instead of
  appending to the input.
- **Free-set early stopping never fired.** The plateau test required a strict
  improvement over the best score, but the free-set CC of a density-modified map
  creeps up by ~0.005 per cycle indefinitely. It is now an absolute-tolerance
  plateau test (`--min-improve`, default 0.02) with arg-max restore.
- **`mask_radius`/`--radius` documentation** claimed the value was not forwarded
  to the engine; it is.

### Changed

- Default cycle budget is `--max-cycles 8` (was 12). True phase error on the
  bundled cases is flat from cycle 4–6 onwards while the FOM keeps rising, so the
  extra cycles cost time and buy nothing.
- Histogram matching now evaluates the ordinal map only on the protein support
  and uses lookup tables instead of two `np.interp` calls over the whole grid:
  ~37% faster per cycle with identical results (checked to 0.02°).
- `_build_modifier` (mask + histogram LUT) moved into the shared engine; v1 and
  v2 no longer carry two copies that can drift apart.

### Added

- **Bundled histogram library** (`oasis_dm/data/protein_hist.lib`), generated
  from public PDB entries (1crn, 1ubq, 2vb1, 1ake) by `tools/make_hist_lib.py`,
  so histogram matching works without a CCP4 installation. CCP4's own `hist.lib`
  still takes priority when it can be found. Measured on the bundled cases it is
  on par with CCP4's library (mean weighted phase error 49.98° vs 50.63°), better
  on two cases and worse on one.
- `--hist-mode protein` — condition the target distribution on the protein volume
  fraction, consistent with matching the protein region only.
- `--per-shell-d` — use the per-resolution-shell σ_A in the `FWT`/`FDM` output
  instead of a single global value.
- `tools/benchmark.py` — reproducible benchmark against reference phases, with an
  optional CCP4 `cparrot` comparison. Every number in the README comes from it.
- Unit tests for the numerics (von Mises round trip, HL centroid, FFT round trip,
  histogram library), output idempotence and run determinism; GitHub Actions CI
  on Python 3.9–3.13 (Linux, macOS, Windows).

### Removed

- An experimental "calibrated" combination weight that derived the von Mises
  precision from the free-set σ_A of the omit map. It is not a valid estimator —
  an omit map still contains ~95% of the reflections, so the amplitude ratio stays
  near 1 and the resulting FOM saturates at 0.98. Removed rather than shipped.

## [2.0.0] — 2026-11-21

First public release: v1 engine (solvent masking, histogram matching, HL
recombination, γ correction) plus the v2 wrapper (omit-map cross-validation,
free-set early stopping), packaged as `oasis_dm` with the `handdm_v2.py`
compatibility shim.
