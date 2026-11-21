# oasis_dm

In-house density modification (phase improvement) for macromolecular crystallography.

`oasis_dm` is a self-contained implementation of the density-modification step used in the
[SAD / SADR pipeline](https://github.com/WeiDing535/OASIS4.3) — solvent masking, histogram
matching, Hendrickson–Lattman phase recombination and free-set driven early stopping. It calls
no CCP4 program and needs no CCP4 installation.

> **What this is not.** `oasis_dm` is a *density modification* engine, not a phasing engine and
> not a drop-in replacement for CCP4 `cparrot`/`parrot`. It improves phases you already have.
> See [Known limitations](#known-limitations) before wiring it into a pipeline.

---

## Install

```bash
pip install .
# or, for development
pip install -e ".[test]"
```

Requires Python ≥ 3.9 and `gemmi`, `numpy`, `scipy`.

### The histogram library

Histogram matching uses `hist.lib`, a **CCP4 data file that is not redistributed with this
package**. Point the program at your own copy in any of these ways:

```bash
oasis-dm --hist-lib /path/to/hist.lib ...   # explicit
export OASIS_DM_HIST_LIB=/path/to/hist.lib  # environment variable
export CLIBD=/path/to/ccp4/lib/data        # CCP4's own lib/data directory
```

If none is found, the program says so clearly and you can run without histogram matching
using `--no-hist` (results are worse but the run completes).

---

## Usage

```bash
oasis-dm --mtzin input.mtz --mtzout output.mtz \
         --fo FOBS --sigf SIGFOBS --phib PHIB --fom FOM --free FreeR \
         --solc 0.5 --cycles 5 \
         --stats results/stats.json
```

`python -m oasis_dm ...` is equivalent. A quick run on a bundled case, scored against
reference phases:

```bash
./examples/run_example.sh 7e1d
```

### Input columns

The input MTZ must contain observed amplitudes and a starting phase estimate. Column names are
set on the command line and default to the CCP4 conventions:

| Flag | Default | Meaning |
| --- | --- | --- |
| `--fo` | `FOBS` | observed amplitude |
| `--sigf` | `SIGFOBS` | standard deviation of `FOBS` |
| `--phib` | `PHIB` | starting phase, degrees |
| `--fom` | `FOM` | starting figure of merit for `PHIB` |
| `--hl` | `HLA,HLB,HLC,HLD` | Hendrickson–Lattman coefficients (optional) |
| `--free` | `FreeR` | free-set flag column; `--free-marker` selects the free value (default `0`) |
| `--truth` | `PHIC` | reference phase, used only for self-assessment |

### Output columns

The output MTZ carries the input columns plus:

| Column | Type | Meaning |
| --- | --- | --- |
| `DM_FDM` | F | density-modification map coefficients |
| `DM_PHIDM` | P | improved phase, degrees |
| `DM_FOMDM` | W | figure of merit of `DM_PHIDM` |
| `DM_FWT` | F | Fourier map coefficients from the improved phases |
| `DM_PHWT` | P | phases of `DM_FWT` |

`--stats` writes a JSON record of the run (per-cycle figures of merit, free-set statistics,
SigmaQ estimates and the settings used).

---

## Scoring a run

`tools/evaluate.py` compares a result against reference phases and reports phase-error
statistics (weighted mean phase error, mean `cos Δφ`, a σ_A-like correlation coefficient),
map correlation coefficients, and per-resolution-shell breakdowns. Work and free sets are
reported separately.

```bash
python tools/evaluate.py \
  --case tests/data/7e1d.mtz --truth-col PHIC --fo-col FOBS --fom-col FOM --free-col FreeR \
  --methods "initial=tests/data/7e1d.mtz:PHIB:FOM,v2=out/7e1d_v2.mtz:DM_PHIDM:DM_FOMDM" \
  --out out/7e1d_eval
```

`python tools/evaluate.py --selftest` runs a synthetic self-check.

---

## Verified results

All numbers below were produced with the code in this repository.

### Phase accuracy against reference phases

Bundled case `7e1d` (P 21 21 21, 2.50 Å, 6447 reflections, 340 free):

| Metric | Input (`PHIB`/`FOM`) | `oasis_dm` output |
| --- | ---: | ---: |
| weighted mean phase error (°) | 56.08 | **45.96** |
| mean phase error (°) | 75.00 | **61.63** |
| mean cos Δφ | 0.460 | **0.584** |
| σ_A-like C | 0.541 | **0.720** |
| map CC (whole cell) | 0.406 | **0.602** |

### End-to-end model building

Measured in the OASIS4.3 pipeline on the `6zpp` calibration case (P 1 21 1, 2.50 Å), all with
the same Buccaneer/REFMAC recipe:

| Configuration | R-work | R-free |
| --- | ---: | ---: |
| CCP4 Parrot (baseline) | 0.4230 | 0.4392 |
| `oasis_dm` + map conditioning | 0.3960 | 0.4207 |
| `oasis_dm` + refinement bridge | **0.3347** | **0.3485** |

These end-to-end numbers **depend on downstream CCP4 programs** (REFMAC and Buccaneer). They
measure `oasis_dm` as the density-modification stage inside a working pipeline, not as a
standalone builder.

---

## Known limitations

Read this before using `oasis_dm` to replace an existing density-modification step.

1. **The raw `oasis_dm` map is often not directly traceable by Buccaneer.** On `6zpp`, the
   phases improve but the initial build from the raw DM map traces few residues. Two workarounds
   are known to work, and they are alternatives, not a stack:
   - *map conditioning* — rescale `FWT` and mix `PHWT` toward the probability phase before
     building (the pipeline's `--oasis-dm-map-scale` / `--oasis-dm-map-phase-mix` options);
   - *refinement bridge* — refine a partial model against the frozen DM phases and rebuild from
     the resulting σ_A model-phase map (the pipeline's `--oasis-dm-rebuild-cycles` option).
2. **A single calibration case is not a benchmark.** The published numbers are one case for
   phase accuracy and one case for end-to-end building. On `7e1d`, a lower combined weight
   (`--cmbwt 0.1`) improves Buccaneer building but degrades the full pipeline; the default
   `--cmbwt 0.4` is kept for that reason.
3. **Parameter sensitivity is not fully mapped.** `mask_radius` is currently not forwarded to
   the engine, and anisotropic correction is not implemented.
4. **Free-set handling.** The engine uses the free-set flag to drive early stopping and, with
   `--omit-free`, for unbiased SigmaQ estimation. It does not use free reflections as
   additional phase observations.

---

## Repository layout

```
src/oasis_dm/        the engine
  core.py            v1 engine and shared numerics (masking, HL recombination, FOM handling)
  dm.py              v2 engine: omit-based unbiased SigmaQ + free-set early stopping
  cli.py             console entry point
  paths.py           runtime data lookup (hist.lib)
tools/evaluate.py    scoring against reference phases
tests/               smoke tests and small bundled cases
examples/            runnable example
handdm_v2.py         compatibility shim for pipelines that resolve the engine by file name
```

### Compatibility with existing deployments

Earlier OASIS releases resolved the in-house engine as a script named `handdm_v2.py`. That name
is kept at the repository root so an existing pipeline keeps working unchanged; it simply
forwards to the packaged CLI.

---

## Data and licensing

- **Code**: MIT — see [LICENSE](LICENSE).
- **Bundled test cases** (`tests/data/*.mtz`): derived from public PDB entries (`3o2e`,
  `6zpp`, `7e1d`) with simulated SAD observations. They contain `FOBS`, `SIGFOBS`, `FreeR`,
  `PHIB`, `FOM`, `PHIC` and the Hendrickson–Lattman coefficients `HLA`–`HLD`.
- **`hist.lib`** is a CCP4 data file and is **not** included. Obtain it from your own CCP4
  installation.

## Citation

If you use `oasis_dm` in published work, please cite the OASIS4.3 pipeline paper and state
that density modification was performed with `oasis_dm` (DOI to be assigned).
