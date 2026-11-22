# oasis_dm

**Density modification (phase improvement) for macromolecular crystallography, without CCP4.**

[![tests](https://github.com/WeiDing535/oasis-dm/actions/workflows/tests.yml/badge.svg)](https://github.com/WeiDing535/oasis-dm/actions/workflows/tests.yml)
[![Python](https://img.shields.io/badge/python-3.9%20%7C%203.10%20%7C%203.11%20%7C%203.12%20%7C%203.13-blue)](pyproject.toml)
[![License: Academic Use](https://img.shields.io/badge/license-academic%20use-blue)](LICENSE)

`oasis_dm` is a self-contained implementation of the density-modification step used in the
[SAD / SADR pipeline](https://github.com/WeiDing535/OASIS4.3): solvent masking, histogram
matching, Hendrickson–Lattman phase recombination, perturbation-γ bias correction and
free-set driven stopping. It calls no CCP4 program and needs no CCP4 installation.

On three bundled SAD cases it lands within a couple of degrees of CCP4 `cparrot` — see
[Benchmark](#benchmark). It is **not** a drop-in replacement: read
[Known limitations](#known-limitations) before wiring it into a pipeline.

---

## Install

```bash
pip install .            # or: pip install -e ".[test]" for development
```

Requires Python ≥ 3.9 and `gemmi`, `numpy`, `scipy`. A histogram library is bundled, so
nothing else is needed to run.

## Quick start

```bash
oasis-dm --mtzin input.mtz --mtzout output.mtz \
         --fo FOBS --sigf SIGFOBS --phib PHIB --fom FOM --free FreeR \
         --solc 0.5

# bundled case, scored against reference phases
./examples/run_example.sh 7e1d
```

`python -m oasis_dm ...` is equivalent. The engine stops on its own when the free-set score
plateaus; `--no-auto-stop --cycles N` forces a fixed number of cycles.

For a multi-copy ASU with known non-crystallographic symmetry, add the operator file:

```bash
oasis-dm --mtzin input.mtz --mtzout output.mtz ... \
         --ncs ncs_operators.txt     # or --ncs-auto to detect + gate it
```

## Desktop app

A native window (PySide6 / Qt 6) that maps one-to-one onto the command line: pick the
input MTZ, it auto-detects the Fo / σ(Fo) / PHIB / FOM / free columns and suggests an
output name; set solvent content, cycles, early stop, mask, histogram target, phase
extension, NCS mode and the advanced SigmaQ options; the exact `oasis-dm` command line
is shown live before you press Run. The engine then runs as a subprocess and its log
streams into the window, with free-set cc plotted per cycle. Stop kills it cleanly.

```bash
pip install ".[gui]"     # once: adds PySide6-Essentials (LGPL)
oasis-dm-gui             # or: python -m oasis_dm.gui
```

**Double-clickable macOS app.** On macOS you can build a self-contained bundle that
needs no Python installation at all — it embeds its own venv with PySide6:

```bash
python3 tools/make_macos_app.py     # -> dist/OasisDM.app
open dist/OasisDM.app               # or double-click it in Finder
```

Why PySide6 and not Tkinter? The Tcl/Tk that ships with many system Pythons is broken —
on the reference machine (macOS 15, Xcode toolchain Python 3.9.6) merely calling
`tk.Tk()` aborts the process. PySide6 installs in one `pip` command, works there, and is
LGPL-licensed. PySide6 is an optional extra: the package and the CLI never need it, and
the test suite skips the window tests when it is absent.

---

## NCS averaging

`--ncs FILE` averages the density over NCS-related copies inside every cycle, after solvent
flattening and histogram matching. The file lists one operator per three lines, fractional
coordinates, `x' = R_f x + t_f` (the format `tools/make_test_case.py --ncs-out`
produces):

```
# A -> B
r11 r12 r13 t1
r21 r22 r23 t2
r31 r32 r33 t3
```

Improper rotations (det = −1) are rejected; proper rotations of any order are accepted, and
open operators (A→B without B→A) are closed by adding the inverse walk.

What the averager actually does, and why:

* **Lattice variants.** An operator fitted on unwrapped coordinates obeys
  `O(f + L) = O(f) + R·L`, so the 27 unit-cell offsets `t + R·L` are all valid placements of
  the same operator. Each mask region uses the variant that actually lands on density; the
  grid cannot be pre-rolled because `R·L` is not generally an integer offset.
* **Voronoi partition.** Voxels whose local map correlation clearly favours one variant
  (double box-smoothed cc > 0.30) act as seeds; every other mask voxel adopts the variant of
  the nearest seed. This follows the molecular copies even when the protein mask is one
  connected blob.
* **Confidence weights.** The blend weight is the smoothed local correlation of the assigned
  variant: ≤ 0.10 → no averaging, ≥ 0.25 → full replacement. Regions where no partner can be
  found keep the un-averaged density, so weak or wrong operators degrade the map gracefully
  instead of smearing it.
* **Space-group closure.** The operator set is closed under crystallographic conjugation, so
  every crystallographic image of a copy is averaged with its own mate and the result stays
  space-group symmetric.

`--ncs-auto` leaves averaging off for the first cycles, then detects an
operator from the current map (a rotation scan on a tight core mask, Nelder–Mead refinement
of rotation + translation) and enables it **only if the free-set correlation beats plain
density modification by a clear margin** (+0.03; the free metric is amplitude-only and
cannot see phase damage, so ties fall back to plain); otherwise the run continues as plain
density modification. It is a conservative safety net, not a substitute for the operator
file: on the bundled cases it declines to average every time — a spurious operator on the
non-NCS control (`7e1d`), a harmful one on `3i0t_hard`, and on `7tj1` a single
half-detected operator (one of the three true ones, local cc 0.32) whose averaging looked
+0.014 better on the free amplitudes but cost 0.07 cos of truth. The measured `7tj1` win
needs the complete expert operator file.

Measured on the bundled NCS cases (mean cos phase error vs reference, higher is better):

| case | plain | `--ncs` | copies | operator CA rmsd | starting wMPE |
| --- | ---: | ---: | ---: | ---: | ---: |
| 7tj1 | 0.538 | **0.552** | 4 | 0.16–0.46 Å | 58° |
| 3i0t | 0.728 | 0.679 | 2 | 1.02 Å | 53° |
| 3i0t_hard | 0.358 | 0.330 | 2 | 1.02 Å | 70° |

Read the table as the honest caveat: **NCS averaging pays when the operators are
tight and the copies are many**. On 7tj1 (four copies, sub-Å operators) averaging
supplies information that solvent flattening alone cannot and wins +0.014 cos. On
3i0t the two chains differ by 1.02 Å CA rmsd — a loose operator — and forcing the
map onto the NCS average costs more than it gains; the poorer the starting phases,
the deeper the loss (3i0t_hard). That is why `--ncs` is documented as an expert
option and `--ncs-auto` gates on the free set: on the bundled non-NCS control
(7e1d) the gate detects a spurious operator and never enables averaging.

---

## Benchmark

Reproduce everything below with:

```bash
python tools/benchmark.py                   # oasis_dm only, bundled histogram library
python tools/benchmark.py --parrot          # add CCP4 cparrot (needs CCP4 on PATH)
```

Three bundled SAD cases at 2.5 Å, all with 50% solvent, all starting from the same phase
probability. `cparrot` is run with `-solvent-flatten -histogram-match -cycles 5
-solvent-content 0.5`, its HL coefficients converted to a centroid phase and FOM, and scored
exactly like `oasis_dm`'s `DM_PHIDM`/`DM_FOMDM`.

Weighted mean phase error (°) — lower is better:

| case | input | `oasis_dm` (bundled lib) | `oasis_dm` (CCP4 `hist.lib`) | CCP4 `cparrot` |
| --- | ---: | ---: | ---: | ---: |
| 7e1d (P 2₁2₁2₁, 6447 refl) | 56.08 | **44.35** | 46.33 | 47.91 |
| 6zpp (P 1 2₁ 1, 5494 refl) | 59.55 | 57.96 | 56.55 | **55.83** |
| 3o2e (P 4₁2₁2, 2805 refl) | 58.85 | **47.64** | 49.00 | 52.37 |
| mean | 58.16 | **49.98** | 50.63 | 52.04 |

The two `oasis_dm` columns are the same code with different histogram targets: the bundled
library (what you get on a machine without CCP4) and CCP4's own `hist.lib`. They differ by
about a degree in each direction — neither dominates — which is why the target library is
worth choosing deliberately.

Map correlation coefficient against the reference map (bundled library / CCP4 lib / cparrot):

| case | map CC |
| --- | --- |
| 7e1d | 0.406 → **0.629** / 0.595 / 0.569 |
| 6zpp | 0.298 → 0.385 / 0.393 / **0.423** |
| 3o2e | 0.437 → **0.608** / 0.586 / 0.557 |

Three honest observations from that table:

* **Phase accuracy is competitive**: `oasis_dm` is better on two of three cases and slightly
  worse on one. It is not uniformly better, and three cases are not a benchmark suite.
* **The reported FOM is optimistic.** With the bundled library it runs 0.08–0.29 above the true
  ⟨cos Δφ⟩ (0.681 vs 0.605, 0.721 vs 0.434, 0.680 vs 0.550); `cparrot` runs 0.01–0.22 above.
  `DM_FOMDM` is usable as a *relative* weight between reflections, not as an absolute accuracy
  estimate — see [Known limitations](#known-limitations).
* **The bundled histogram library is a substitute.** It happens to beat CCP4's on two of these
  three cases and lose on the third; do not read that as a general ranking.

Ablations on the same cases (weighted mean phase error °, bundled library):

| configuration | 7e1d | 6zpp | 3o2e |
| --- | ---: | ---: | ---: |
| full engine | 44.35 | 57.96 | 47.64 |
| without histogram matching (`--no-hist`) | 51.85 | 57.04 | 57.48 |
| without solvent flattening (`--no-solv`) | 56.06 | 62.22 | 60.54 |

Histogram matching and solvent flattening are both essential; histogram matching is the single
largest contributor, which is why the target distribution matters (next section).

---

## The histogram library

Histogram matching needs a resolution-dependent target distribution of protein density. CCP4
ships one as `$CLIBD/hist.lib`; that file is **not** redistributable, so this package looks
for a target in this order:

1. `--hist-lib PATH`
2. `$OASIS_DM_HIST_LIB`
3. `$CLIBD/hist.lib` (CCP4's own library — preferred if you have it)
4. the bundled `oasis_dm/data/protein_hist.lib`

The bundled file is generated from public PDB entries by `tools/make_hist_lib.py`: each entry
is reduced to amino-acid residues, structure factors are computed with real atomic scattering
factors, truncated at each target resolution and FFTed, and the protein-volume histogram is
averaged over entries. It is a substitute, not a copy of CCP4's curated database.

```bash
python tools/make_hist_lib.py --out src/oasis_dm/data/protein_hist.lib   # regenerate
oasis-dm ... --no-hist        # skip histogram matching entirely
```

---

## Inputs and outputs

The input MTZ must contain observed amplitudes and a starting phase estimate. Column names are
command-line options defaulting to CCP4 conventions:

| Flag | Default | Meaning |
| --- | --- | --- |
| `--fo` | `FOBS` | observed amplitude |
| `--sigf` | `SIGFOBS` | standard deviation of `FOBS` |
| `--phib` | `PHIB` | starting phase, degrees (radians are detected and converted) |
| `--fom` | `FOM` | starting figure of merit |
| `--hl` | `HLA,HLB,HLC,HLD` | HL coefficients; used instead of `PHIB`/`FOM` when present |
| `--free` | `FreeR` | free-set flag column; `--free-marker` selects the free value (default `0`) |
| `--truth` | `PHIC` | reference phase, used only for self-assessment |

Output columns (input columns are preserved):

| Column | Type | Meaning |
| --- | --- | --- |
| `DM_FDM` | F | density-modification map coefficients |
| `DM_PHIDM` | P | improved phase, degrees |
| `DM_FOMDM` | W | figure of merit of `DM_PHIDM` (optimistic, see limitations) |
| `DM_HLADM…HLDDM` | A | HL coefficients of the combined phase probability |
| `DM_FWT` / `DM_PHWT` | F / P | σ_A-weighted map coefficients (what a builder should read) |
| `DM_DELFWT` / `DM_PHDELFWT` | F / P | difference map coefficients |
| `DM_FCDM` / `DM_PHICDM` | F / P | modified-map structure factors |

`--stats FILE` writes a JSON record of the run (per-cycle figures of merit, free-set
statistics, σ_A and Σ_Q per shell, and every setting used).

---

## How it works

Each cycle:

1. **Best map** — centroid phase and FOM of the current HL distribution (numerical quadrature,
   so bimodal SAD priors are handled), projected onto symmetry-allowed values for centric
   reflections; FFT to a P1 grid.
2. **Solvent envelope** — Abrahams–Leslie local-variance mask (`--mask variance`, default) or
   Wang rolling sphere (`--mask wang`), thresholded at `--solc` with a soft sigmoid edge.
3. **Histogram matching** — ordinal (rank) mapping of the protein-region density onto the
   target distribution, with Parrot's mild ordinal squash.
4. **Perturbation γ** — the nonlinearity of the modification is measured by perturbing the map
   and taking the regression slope; the bias it introduces is divided out.
5. **Recombination** — the modified map contributes a von Mises term with weight
   `2·Fo·D·Fc/Σ_Q` (Sim/σ_A), damped by `--cmbwt` and added to the *original* HL coefficients
   (`--anchor exp`), so information cannot compound across cycles.
6. **Cross-validation** — a second map is built from the same mask and lookup table with the
   free reflections deleted; the free-set amplitude CC drives stopping, and the best cycle is
   restored.

## Knobs worth touching

| Option | Default | When to change it |
| --- | --- | --- |
| `--solc` | 0.5 | set it from the Matthews coefficient; ±0.05 changes little, ±0.15 hurts |
| `--mask` / `--radius` | variance, 4 Å | `wang`/8 Å behaves like `dm`'s SOLMASK |
| `--cmbwt` | 0.4 | lower (0.2–0.3) trusts the modified map less |
| `--max-cycles` | 8 | true phase error is flat from ~cycle 5 |
| `--hist-mode` | `cell` | `protein` conditions the target on `--solc` |
| `--per-shell-d` | off | per-shell σ_A in `FWT` (no measurable effect on the bundled cases) |
| `--ncs FILE` | off | NCS operator file; averaging runs from cycle 1 — know your operators |
| `--ncs-auto` | off | detect an operator at cycle ≥ 2 and enable it only if the free-set score improves |

Everything else — `--x-max`, `--n-shells`, `--radius-decay`, `--sample-rate` — is saturated on
the bundled cases: sweeping them moves the weighted phase error by less than 1°.

---

## Known limitations

1. **The reported FOM is optimistic** by 0.08–0.29 (see the benchmark table). The Sim/σ_A weight
   assumes the only error in the modified map is unmodelled structure, which is false when the
   input phases are 40–60° off. `--cmbwt 0.4` damps this; it does not remove it. An attempt to
   calibrate the weight from the free-set σ_A of the omit map does **not** work — the omit map
   still contains ~95% of the reflections, so the amplitude ratio stays near 1 — and was
   removed rather than shipped (see `CHANGELOG.md`).
2. **Three cases, one resolution.** All bundled cases are 2.5 Å SAD data with ~50% solvent.
   Nothing here is calibrated for low-resolution SAD, SIR/MIR, or very high solvent content.
3. **NCS averaging is operator-limited.** The gain is set by the quality of the
   operators, not by how poor the phases are: tight sub-Å operators on many
   copies win (7tj1), a loose operator loses, and the poorer the starting
   phases the deeper the loss. The correlation-weighted blend protects the map
   (regions where no partner is found keep un-averaged density) but cannot
   rescue a wrong or loose operator. Only proper rotations are supported, and
   operators stay fixed during the run (they are not re-refined against the
   improving map).
4. **Free-set handling.** The free flag drives stopping and the omit map; free reflections are
   not used as extra phase observations. Below ~500 free reflections the omit-based σ_Q path
   (`--sigmaq-source omit`) falls back to a work-dominated estimate.
5. **`--anchor cumulative` is available but not recommended** — it inflates the FOM to ~0.88
   and does not improve phases.
6. **The bundled histogram library is a substitute** for CCP4's curated one. If you have CCP4,
   let it find `$CLIBD/hist.lib`.

---

## Repository layout

```
src/oasis_dm/
  core.py            v1 engine + shared numerics (mask, histogram LUT, HL recombination, I/O)
  dm.py              v2 engine: omit map, free-set stopping, phase extension, CLI
  ncs.py             NCS averaging (operator I/O, lattice variants, auto-detection)
  cli.py             console entry point
  gui.py             desktop GUI (PySide6 window, engine subprocess)
  paths.py           histogram-library resolution
  data/protein_hist.lib   bundled target distribution
tools/
  benchmark.py       benchmark vs reference phases (+ optional cparrot)
  evaluate.py        scoring against reference phases (plots optional)
  make_hist_lib.py   regenerate the bundled histogram library
  make_test_case.py  build a bundled-style test case (+ NCS operator files) from a PDB entry
  make_macos_app.py  build a double-clickable dist/OasisDM.app with an embedded venv
  run_tests.py       test runner that needs no pytest
tests/               smoke tests, numerics tests, bundled cases (incl. NCS cases)
examples/            runnable example
handdm_v2.py         compatibility shim for pipelines that resolve the engine by file name
```

### Compatibility with existing deployments

Earlier OASIS releases resolved the in-house engine as a script named `handdm_v2.py`. That name
is kept at the repository root so an existing pipeline keeps working unchanged; it forwards to
the packaged CLI.

---

## Data and licensing

- **Code** is free for noncommercial academic research and teaching under the [Academic Use License](LICENSE). Commercial use requires a paid license. Contact dingwei@iphy.ac.cn. Earlier MIT releases retain their [MIT terms](LICENSE-MIT-PREVIOUS.txt).
- **Bundled test cases** (`tests/data/*.mtz`): derived from public PDB entries (`3o2e`,
  `6zpp`, `7e1d`) with simulated SAD observations. They contain `FOBS`, `SIGFOBS`, `FreeR`,
  `PHIB`, `FOM`, `PHIC` and HL coefficients `HLA`–`HLD`.
  *Note*: the `6zpp` truth column was corrected in 2.1.0 — the original was offset by 180°·k
  along **b** relative to `PHIB`, which made the case look catastrophic. See `CHANGELOG.md`.
- **`oasis_dm/data/protein_hist.lib`**: generated here from public PDB entries, MIT like the
  rest of the package.
- **CCP4 `hist.lib`** is not included; obtain it from your own CCP4 installation.

## Citation

If you use `oasis_dm` in published work, cite the OASIS4.3 pipeline paper and state that
density modification was performed with `oasis_dm` (see [CITATION.cff](CITATION.cff)).

## Development

```bash
pip install -e ".[test]"
python tools/run_tests.py        # no pytest needed
python -m pytest tests -q        # if pytest is available
python tools/benchmark.py        # benchmark, ~30 s on the bundled cases
```
