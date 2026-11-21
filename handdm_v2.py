#!/usr/bin/env python3
"""Compatibility entry point for existing deployments that call `handdm_v2.py`.

Earlier OASIS releases resolved the in-house density-modification engine as a
script named ``handdm_v2.py``. That name is kept here so an existing pipeline
keeps working unchanged; it simply forwards to the packaged CLI.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from oasis_dm.cli import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
