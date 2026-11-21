#!/usr/bin/env python3
"""Run the smoke tests without needing pytest installed.

    python tools/run_tests.py

`pytest tests/` works too when pytest is available; this runner exists so the
package can be validated in a minimal or offline environment.
"""
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))


def main():
    import test_cli

    names = sorted(n for n in dir(test_cli) if n.startswith("test_"))
    failures = []
    for name in names:
        fn = getattr(test_cli, name)
        if not callable(fn):
            continue
        try:
            fn()
            print("PASS %s" % name)
        except Exception:  # noqa: BLE001
            failures.append(name)
            print("FAIL %s" % name)
            traceback.print_exc()
    print("\n%d/%d passed" % (len(names) - len(failures), len(names)))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
