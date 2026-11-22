#!/usr/bin/env python3
"""Build a double-clickable macOS app bundle for the oasis-dm GUI.

Usage (from the repository root)::

    python3 tools/make_macos_app.py

The result is ``dist/OasisDM.app``.  Double-clicking it (or ``open
dist/OasisDM.app``) starts the desktop GUI — no Terminal window, no
PYTHONPATH, nothing to install first.  The bundle embeds its own virtual
environment with PySide6 and this package installed, so it keeps working
independently of whatever happens to the system Python.

The bundle is unsigned: fine for yourself and for colleagues copying it
within a lab (Apple only quarantines files downloaded from the internet,
not apps built locally).  For wider distribution, sign it with
``codesign`` and consider notarisation.

Implementation notes:

- The venv is created from a full Python (default: the interpreter running
  this script, or ``--python``).  A venv created from a real framework
  build relocates fine, because ``bin/python`` finds ``pyvenv.cfg`` next
  to itself at runtime; only the base-Python path inside ``pyvenv.cfg``
  is absolute, and that one is stable system-wide.
- We deliberately do NOT use the system Tcl/Tk anywhere: on several macOS
  + Xcode-Python combinations ``tk.Tk()`` aborts the process outright.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
APP = "OasisDM"

PLIST = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN"
  "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
  <key>CFBundleName</key>              <string>{app}</string>
  <key>CFBundleDisplayName</key>       <string>oasis-dm</string>
  <key>CFBundleIdentifier</key>        <string>com.weiding.oasis-dm</string>
  <key>CFBundleVersion</key>           <string>{version}</string>
  <key>CFBundleShortVersionString</key><string>{version}</string>
  <key>CFBundlePackageType</key>       <string>APPL</string>
  <key>CFBundleExecutable</key>        <string>{app}</string>
  <key>LSMinimumSystemVersion</key>    <string>11.0</string>
  <key>NSHighResolutionCapable</key>   <true/>
  <key>NSSupportsAutomaticTermination</key><false/>
  <key>NSSupportsSuddenTermination</key>   <false/>
</dict>
</plist>
"""

LAUNCHER = """#!/bin/bash
# oasis-dm GUI launcher: runs the venv python embedded in this bundle.
DIR="$(cd "$(dirname "$0")/.." && pwd)"
exec "$DIR/Resources/venv/bin/python" -m oasis_dm.gui "$@"
"""


def run(cmd, **kw):
    print("+", " ".join(str(c) for c in cmd))
    subprocess.run([str(c) for c in cmd], check=True, **kw)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--python", default=sys.executable,
                    help="base python for the embedded venv "
                         "(default: the one running this script)")
    ap.add_argument("--dest", default="dist",
                    help="output directory (default: dist/)")
    args = ap.parse_args(argv)

    from oasis_dm import __version__          # needs PYTHONPATH=src

    dest = (REPO / args.dest / (APP + ".app")).resolve()
    if dest.exists():
        shutil.rmtree(dest)
    macos = dest / "Contents" / "MacOS"
    resources = dest / "Contents" / "Resources"
    macos.mkdir(parents=True)
    resources.mkdir(parents=True)

    exe = macos / APP
    exe.write_text(LAUNCHER)
    exe.chmod(0o755)
    (dest / "Contents" / "Info.plist").write_text(
        PLIST.format(app=APP, version=__version__))

    print("building", dest)
    run([args.python, "-m", "venv", resources / "venv"])
    vpy = resources / "venv" / "bin" / "python"
    run([vpy, "-m", "pip", "install", "--quiet", "--upgrade", "pip"])
    run([vpy, "-m", "pip", "install", "--quiet",
         "PySide6-Essentials>=6.5"])
    run([vpy, "-m", "pip", "install", "--quiet", REPO])
    # import check: gui module must be present and Qt loadable
    run([vpy, "-c",
         "import PySide6.QtWidgets, oasis_dm, oasis_dm.gui; "
         "print('bundle check ok: oasis-dm', oasis_dm.__version__)"])
    print("\ndone: %s   (double-click it, or: open %s)" % (dest, dest))
    return 0


if __name__ == "__main__":
    sys.exit(main())
