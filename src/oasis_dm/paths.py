"""Locate runtime data without hard-coding a platform path.

`hist.lib` is a CCP4 data file and is deliberately not redistributed with this
package. Point `--hist-lib` at your own copy, set ``OASIS_DM_HIST_LIB``, or let
the resolver find it through ``$CLIBD`` (the CCP4 ``lib/data`` directory).
"""
from __future__ import annotations

import os
from pathlib import Path

ENV_VAR = "OASIS_DM_HIST_LIB"
FILENAME = "hist.lib"

MISSING = (
    "no histogram library found: pass --hist-lib PATH, set "
    + ENV_VAR
    + ", or set CLIBD to your CCP4 lib/data directory. hist.lib is a CCP4 data "
      "file and is not redistributed with this package; pass --no-hist to run "
      "without histogram matching."
)


def resolve_hist_lib(explicit: str = "") -> str:
    """Return a readable histogram-library path or raise a clear error."""
    if explicit:
        path = Path(explicit).expanduser()
        if path.is_file():
            return str(path)
        raise ValueError("--hist-lib is not a file: " + str(path))
    candidates = []
    env = os.environ.get(ENV_VAR)
    if env:
        candidates.append(Path(env).expanduser())
    clibd = os.environ.get("CLIBD")
    if clibd:
        candidates.append(Path(clibd) / FILENAME)
    here = Path(__file__).resolve()
    for base in (Path.cwd(), here.parent, here.parent.parent.parent):
        candidates.append(base / FILENAME)
    for candidate in candidates:
        if candidate.is_file():
            return str(candidate)
    raise ValueError(MISSING)
