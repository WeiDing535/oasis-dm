"""Console entry point for ``oasis-dm`` and ``python -m oasis_dm``."""
from __future__ import annotations

import sys

from .dm import main as _run


def main(argv=None) -> int:
    return _run(argv)


if __name__ == "__main__":
    sys.exit(main())
