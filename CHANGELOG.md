# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[semantic versioning](https://semver.org/).

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
