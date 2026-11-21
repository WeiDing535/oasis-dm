"""Locate the histogram-matching target distribution.

Priority order:

1. ``--hist-lib PATH`` (explicit);
2. ``$OASIS_DM_HIST_LIB``;
3. ``$CLIBD/hist.lib`` -- CCP4's own library, the reference target;
4. the library bundled with this package (``oasis_dm/data/protein_hist.lib``),
   derived from public PDB entries by ``tools/make_hist_lib.py``;
5. ``hist.lib`` in the working directory or next to the source tree.

CCP4's ``hist.lib`` is not redistributed here, so the bundled file is a
self-generated substitute: same format and conventions, but averaged over a
handful of structures instead of a curated database.  ``oasis_dm`` logs which one
it used.
"""
from __future__ import annotations

import os
from pathlib import Path

ENV_VAR = "OASIS_DM_HIST_LIB"
FILENAME = "hist.lib"
BUNDLED = "protein_hist.lib"

MISSING = (
    "no histogram library found: pass --hist-lib PATH to a hist.lib file, or set "
    f"{ENV_VAR} (or CLIBD to your CCP4 lib/data directory). The bundled fallback "
    "(hist.lib derived from public PDB entries) is missing from this "
    "installation; pass --no-hist to run without histogram matching."
)


def bundled_path() -> Path:
    return Path(__file__).resolve().parent / "data" / BUNDLED


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
    candidates.append(bundled_path())
    here = Path(__file__).resolve()
    for base in (Path.cwd(), here.parent, here.parent.parent.parent):
        candidates.append(base / FILENAME)
    for candidate in candidates:
        if candidate.is_file():
            return str(candidate)
    raise ValueError(MISSING)
