# Bundled test cases

Three small SAD cases used by the smoke tests, the example script and the
benchmark. Each is derived from a public PDB entry with simulated SAD
observations and reference phases, and carries the same columns:

`H K L FOBS SIGFOBS FreeR PHIB FOM PHIC HLA HLB HLC HLD`

| file | PDB entry | space group | d_min | reflections | free set | starting wMPE |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| `7e1d.mtz` | 7E1D | P 2₁2₁2₁ | 2.50 Å | 6447 | 340 | 56.1° |
| `6zpp.mtz` | 6ZPP | P 1 2₁ 1 | 2.50 Å | 5494 | 292 | 59.5° |
| `3o2e.mtz` | 3O2E | P 4₁2₁2 | 2.50 Å | 2805 | 143 | 58.8° |

`PHIB`/`FOM` are the starting phase estimate, `HLA`–`HLD` its Hendrickson–Lattman
form, and `PHIC` the reference phase used only for scoring. `FreeR` marks the free
set with the value `0`.

## NCS averaging cases

Two additional 2-fold / 4-copy cases for the NCS tests, same columns and also
2.5 Å simulated SAD data. Each ships with a `_ncs.txt` operator file (Kabsch
fit on Cα atoms, fractional coordinates, `x' = R_f x + t_f`, one operator per
three `r r r t` rows):

| file | PDB entry | space group | d_min | reflections | copies | starting wMPE | NCS |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- |
| `3i0t.mtz` | 3I0T | P 2₁2₁2₁ | 2.50 Å | 19203 | 2 | 53.3° | 179.5° NCS, CA rmsd 1.02 Å |
| `3i0t_hard.mtz` | 3I0T | P 2₁2₁2₁ | 2.50 Å | 19203 | 2 | **70.2°** | same operator, worse starting phases |
| `7tj1.mtz` | 7TJ1 | P 2₁2₁2₁ | 2.50 Å | 18839 | 4 | 58.4° | 3 operators (~180°), CA rmsd 0.16–0.46 Å |

`3i0t_hard.mtz` is the same data as `3i0t.mtz` with a more perturbed starting
phase set — the case where the loose 1.02 Å operator hurts most (the poorer the
phases, the deeper the loss), so `--ncs-auto` must decline to average on it.
`3i0t.mtz` itself is the easy counterpart where averaging is not needed and must
do no harm; `7tj1.mtz` is the positive case (tight operators, four copies) and
`7e1d.mtz` above serves as the non-NCS control for `--ncs-auto`.

Recreate or build new cases (including the operator files) with
`tools/make_test_case.py`.

## 6zpp: corrected truth column (2.1.0)

In the 2.0.0 release, `PHIC` in `6zpp.mtz` was related to `PHIB` by a half-cell
origin shift along **b**, so the reference and the starting phases disagreed
completely (weighted mean phase error 92°, map CC −0.02). That made the case look
as though the engine destroyed the phases. The file now stores
`PHIC + 180·k`, which restores the intended relationship; the starting phases
score 59.5° with map CC 0.298, comparable to the other two cases.

The origin of a P 1 2₁ 1 structure along **b** is not fixed by symmetry, so both
conventions are internally valid — only their *combination* was wrong.

## Scoring

```bash
python tools/evaluate.py --case tests/data/7e1d.mtz \
  --truth-col PHIC --fo-col FOBS --fom-col FOM --free-col FreeR \
  --methods "initial=tests/data/7e1d.mtz:PHIB:FOM,dm=out/7e1d.mtz:DM_PHIDM:DM_FOMDM" \
  --out out/7e1d_eval
```

or, for all cases at once, `python tools/benchmark.py`.
