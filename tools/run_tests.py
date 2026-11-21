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
    import importlib

    modules = sorted(p.stem for p in (ROOT / "tests").glob("test_*.py"))
    names, failures = [], []
    for mod_name in modules:
        mod = importlib.import_module(mod_name)
        for name in sorted(n for n in dir(mod) if n.startswith("test_")):
            fn = getattr(mod, name)
            if not callable(fn):
                continue
            label = "%s.%s" % (mod_name, name)
            names.append(label)
            try:
                fn()
                print("PASS %s" % label)
            except Exception:  # noqa: BLE001
                failures.append(label)
                print("FAIL %s" % label)
                traceback.print_exc()
    print("\n%d/%d passed" % (len(names) - len(failures), len(names)))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
