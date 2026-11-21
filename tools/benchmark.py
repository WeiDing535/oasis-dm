#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Benchmark oasis_dm against reference phases, and optionally against CCP4 cparrot.

Runs the packaged CLI on each bundled case, scores the result with
``tools/evaluate.py`` (weighted mean phase error, mean cos dphi, sigma_A-like C,
map CC) and prints a markdown table.  Every number in the README's benchmark
section comes from this script.

    python tools/benchmark.py                     # bundled cases
    python tools/benchmark.py --cases 7e1d 3o2e
    python tools/benchmark.py --parrot            # also run cparrot, if installed
    python tools/benchmark.py --out out/bench

``--parrot`` needs ``cparrot`` on PATH (or CCP4's ``bin`` on PATH); it is run with
solvent flattening + histogram matching, 5 cycles, 50% solvent -- the closest
equivalent recipe.  cparrot's phase probability is converted to a centroid phase
and FOM before scoring, which is exactly what oasis_dm outputs as
``DM_PHIDM``/``DM_FOMDM``.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CASES = ["7e1d", "6zpp", "3o2e"]
def run(cmd, env=None, cwd=None):
    return subprocess.run(cmd, cwd=str(cwd or ROOT), capture_output=True,
                          text=True, env=env)


def score(case_mtz, method_mtz, phi_col, fom_col, out_prefix, env):
    """Return the metric dict for one method via tools/evaluate.py."""
    r = run([sys.executable, str(ROOT / "tools" / "evaluate.py"),
             "--case", str(case_mtz), "--truth-col", "PHIC",
             "--fo-col", "FOBS", "--fom-col", "FOM", "--free-col", "FreeR",
             "--methods", "m=%s:%s:%s" % (method_mtz, phi_col, fom_col),
             "--out", str(out_prefix)], env=env)
    if r.returncode != 0:
        raise RuntimeError("evaluate.py failed:\n" + r.stderr[-2000:])
    data = json.loads(Path(str(out_prefix) + "_metrics.json").read_text())
    m = data["methods"]["m"]
    g, wf, cc = m["global_"], m["work_free"], m["map_CC"]
    return dict(wMPE=g["wMPE"], MPE=g["MPE"], cos=g["mean_cos_dphi"],
                C=g["sigmaA_like_C"], FOM=g["mean_FOM"],
                CC=cc["whole_cell"], wfree=wf["wMPE_free"])


def centroid_from_hl(mtz_path, env):
    """cparrot's HL columns -> centroid phase/FOM, written to a temp MTZ."""
    code = r'''
import sys, numpy as np, gemmi
sys.path.insert(0, %r)
from oasis_dm.core import hl_centroid
src, out = sys.argv[1], sys.argv[2]
m = gemmi.read_mtz_file(src)
hl = np.column_stack([np.array(m.column_with_label("parrot.ABCD." + c).array, float)
                      for c in "ABCD"])
phi, fom = hl_centroid(hl)
o = gemmi.Mtz(with_base=True)
o.spacegroup, o.cell = m.spacegroup, m.cell
o.add_dataset("parrot")
for c in m.columns[3:]:
    if c.label.startswith("parrot."):
        continue
    o.add_column(c.label, c.type)
o.add_column("PARROT_PHI", "P"); o.add_column("PARROT_FOM", "W")
data = np.column_stack([np.asarray(m.array)[:, [i for i, c in enumerate(m.columns)
                        if not c.label.startswith("parrot.")]],
                        np.degrees(phi), np.clip(fom, 0, 1)]).astype(np.float32)
o.set_data(np.ascontiguousarray(data))
o.write_to_file(out)
''' % str(ROOT / "src")
    tmp = str(Path(os.path.dirname(mtz_path)) / "_parrot_centroid.mtz")
    r = run([sys.executable, "-c", code, str(mtz_path), tmp], env=env)
    if r.returncode != 0:
        raise RuntimeError("centroid conversion failed:\n" + r.stderr[-1500:])
    return tmp, "PARROT_PHI", "PARROT_FOM"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cases", nargs="+", default=DEFAULT_CASES,
                    help="bundled case tags (tests/data/<tag>.mtz)")
    ap.add_argument("--cycles", type=int, default=None,
                    help="fixed cycle count; default uses the engine's auto stop")
    ap.add_argument("--parrot", action="store_true", help="also run CCP4 cparrot")
    ap.add_argument("--out", default="out/bench", help="output directory")
    ap.add_argument("--extra", nargs=argparse.REMAINDER, default=[],
                    help="extra oasis_dm options, e.g. --extra --no-hist "
                         "(must come last)")
    args = ap.parse_args(argv)

    outdir = Path(args.out).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    env = {**os.environ, "PYTHONPATH": str(ROOT / "src")}
    parrot = shutil.which("cparrot") if args.parrot else None
    if args.parrot and not parrot:
        print("note: cparrot not on PATH, skipping the CCP4 comparison")

    rows = {}
    for tag in args.cases:
        case = ROOT / "tests" / "data" / (tag + ".mtz")
        if not case.is_file():
            print("skip %s (missing %s)" % (tag, case))
            continue
        mtz_out = outdir / (tag + "_dm.mtz")
        cmd = [sys.executable, "-m", "oasis_dm", "--mtzin", str(case),
               "--mtzout", str(mtz_out), "--solc", "0.5", "--quiet"]
        if args.cycles:
            cmd += ["--no-auto-stop", "--cycles", str(args.cycles)]
        cmd += [a for a in args.extra if a != "--"]
        r = run(cmd, env=env)
        if r.returncode != 0:
            print("FAILED %s:\n%s" % (tag, r.stderr[-1500:]))
            continue
        rows.setdefault("input", {})[tag] = score(
            case, str(case), "PHIB", "FOM", outdir / (tag + "_in"), env)
        rows.setdefault("oasis_dm", {})[tag] = score(
            case, str(mtz_out), "DM_PHIDM", "DM_FOMDM", outdir / (tag + "_dm"), env)

        if parrot:
            pmtz = outdir / (tag + "_parrot.mtz")
            r = run([parrot, "-mtzin", str(case), "-colin-fo", "FOBS,SIGFOBS",
                     "-colin-phifom", "PHIB,FOM", "-mtzout", str(pmtz),
                     "-solvent-flatten", "-histogram-match", "-cycles", "5",
                     "-solvent-content", "0.5"])
            if r.returncode == 0 and pmtz.is_file():
                tmp, pc, pf = centroid_from_hl(pmtz, env)
                rows.setdefault("cparrot", {})[tag] = score(
                    case, tmp, pc, pf, outdir / (tag + "_parrot"), env)

    if not rows:
        raise SystemExit("nothing ran")

    tags = [t for t in args.cases if any(t in v for v in rows.values())]
    lines = ["| method | " + " | ".join(t for t in tags) + " |",
             "| --- |" + " ---: |" * len(tags)]
    for name in ("input", "oasis_dm", "cparrot"):
        if name not in rows:
            continue
        cells = []
        for t in tags:
            r = rows[name].get(t)
            if r is None:
                cells.append("n/a")
                continue
            cells.append("wMPE %.2f / cos %.3f / CC %.3f / FOM %.3f"
                         % (r["wMPE"], r["cos"], r["CC"], r["FOM"]))
        lines.append("| %s | %s |" % (name, " | ".join(cells)))
    table = "\n".join(lines)
    print()
    print(table)
    (outdir / "benchmark.json").write_text(json.dumps(rows, indent=2))
    (outdir / "benchmark.md").write_text(table + "\n")
    print("\nwrote %s" % (outdir / "benchmark.json"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
