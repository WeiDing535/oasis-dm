"""oasis_dm -- in-house macromolecular density modification (phase improvement).

`oasis_dm` replaces the density-modification step (CCP4 ``dm`` / ``cparrot``)
with a self-contained implementation: solvent masking, histogram matching,
HL phase recombination and free-set diagnostics, with no CCP4 program calls.

Public entry points
-------------------
``oasis-dm`` (console script) or ``python -m oasis_dm`` runs the v2 engine.
``oasis_dm.core`` is the v1 engine and the shared numerical library.
"""

__version__ = "2.0.0"

__all__ = ["__version__"]
