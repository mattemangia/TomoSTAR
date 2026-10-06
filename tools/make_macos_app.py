#!/usr/bin/env python3
# Copyright 2026 Matteo Mangiagalli
# SPDX-License-Identifier: Apache-2.0
"""Builds TomoSTAR.app, a macOS application bundle around the two self-contained executables (Apple silicon and
Intel), and zips it with the execute permissions kept. Double-clicking the application opens a Terminal window in
which `tomostar` is on the PATH and prints its help; Contents/Resources/bin/tomostar runs the executable of the
machine's architecture and can be linked into /usr/local/bin.

Usage (after `dotnet publish -r osx-arm64` and `-r osx-x64` with --self-contained -p:PublishSingleFile=true):
    python3 tools/make_macos_app.py VERSION ARM64_EXE X64_EXE OUT_DIR
Writes OUT_DIR/TomoSTAR.app and OUT_DIR/TomoSTAR-VERSION-macos.zip. Needs Pillow for the icon."""
import os, shutil, sys, time, zipfile
from PIL import Image, ImageDraw

version, arm64, x64, out = sys.argv[1:5]
repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
app = os.path.join(out, "TomoSTAR.app")
shutil.rmtree(app, ignore_errors=True)
C = os.path.join(app, "Contents")
for d in ("MacOS", "Resources/bin/arm64", "Resources/bin/x86_64"):
    os.makedirs(os.path.join(C, d))

def write(path, text, exe=False):
    with open(path, "w", newline="\n") as f:
        f.write(text)
    if exe:
        os.chmod(path, 0o755)

write(os.path.join(C, "Info.plist"), f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
  <key>CFBundleName</key><string>TomoSTAR</string>
  <key>CFBundleDisplayName</key><string>TomoSTAR</string>
  <key>CFBundleIdentifier</key><string>io.github.mattemangia.tomostar</string>
  <key>CFBundleVersion</key><string>{version}</string>
  <key>CFBundleShortVersionString</key><string>{version}</string>
  <key>CFBundlePackageType</key><string>APPL</string>
  <key>CFBundleExecutable</key><string>TomoSTAR</string>
  <key>CFBundleIconFile</key><string>TomoSTAR</string>
  <key>LSMinimumSystemVersion</key><string>14.0</string>
  <key>NSHumanReadableCopyright</key><string>Copyright 2026 Matteo Mangiagalli. Apache License 2.0.</string>
</dict>
</plist>
""")

# The application's executable: open the shell window.
write(os.path.join(C, "MacOS", "TomoSTAR"), """#!/bin/bash
# Opens a Terminal window with tomostar on the PATH.
R="$(cd "$(dirname "$0")/../Resources" && pwd)"
open -a Terminal "$R/TomoSTAR.command"
""", exe=True)

write(os.path.join(C, "Resources", "TomoSTAR.command"), f"""#!/bin/bash
# Terminal session of TomoSTAR.app: tomostar on the PATH, its help, then an interactive shell.
BIN="$(cd "$(dirname "$0")/bin" && pwd)"
export PATH="$BIN:$PATH"
clear
echo "TomoSTAR {version}: seismic travel-time and attenuation tomography from the command line"
echo
tomostar help
echo
echo "tomostar is on the PATH of this window. To use it in every terminal, run once:"
echo "  sudo mkdir -p /usr/local/bin && sudo ln -sf \\"$BIN/tomostar\\" /usr/local/bin/tomostar"
echo
exec "${{SHELL:-/bin/zsh}}" -i
""", exe=True)

# Architecture selector, which also works through a symbolic link.
write(os.path.join(C, "Resources", "bin", "tomostar"), """#!/bin/bash
# Runs the TomoSTAR executable of this Mac's architecture.
p="$0"
while [ -L "$p" ]; do t="$(readlink "$p")"; case "$t" in /*) p="$t" ;; *) p="$(dirname "$p")/$t" ;; esac; done
D="$(cd "$(dirname "$p")" && pwd)"
case "$(uname -m)" in
  arm64) exec "$D/arm64/tomostar" "$@" ;;
  *) exec "$D/x86_64/tomostar" "$@" ;;
esac
""", exe=True)
for src, arch in ((arm64, "arm64"), (x64, "x86_64")):
    dst = os.path.join(C, "Resources", "bin", arch, "tomostar")
    shutil.copy2(src, dst)
    os.chmod(dst, 0o755)
for f in ("LICENSE", "NOTICE", "THIRD_PARTY_NOTICES.md"):
    shutil.copy2(os.path.join(repo, f), os.path.join(C, "Resources", f))

# Icon: rays from a source (star) to receivers (triangles) over an inversion grid.
S = 1024
img = Image.new("RGBA", (S, S), (0, 0, 0, 0))
g = ImageDraw.Draw(img)
g.rounded_rectangle([40, 40, S - 40, S - 40], radius=200, fill=(18, 44, 82, 255))
for k in range(1, 8):
    x = 40 + k * (S - 80) / 8
    g.line([(x, 140), (x, S - 140)], fill=(48, 84, 130, 255), width=4)
    g.line([(140, x), (S - 140, x)], fill=(48, 84, 130, 255), width=4)
src = (512, 790)
recv = [(230, 250), (400, 215), (620, 215), (790, 250)]
for r in recv:
    pts = [(src[0] + (r[0] - src[0]) * t + 60 * t * (1 - t) * (1 if r[0] > src[0] else -1), src[1] + (r[1] - src[1]) * t) for t in [i / 40 for i in range(41)]]
    g.line(pts, fill=(255, 196, 64, 255), width=14, joint="curve")
for r in recv:
    g.polygon([(r[0], r[1] - 46), (r[0] - 42, r[1] + 26), (r[0] + 42, r[1] + 26)], fill=(255, 255, 255, 255))
import math
star = [(src[0] + (60 if i % 2 == 0 else 26) * math.cos(math.pi / 2 + i * math.pi / 5),
         src[1] - (60 if i % 2 == 0 else 26) * math.sin(math.pi / 2 + i * math.pi / 5)) for i in range(10)]
g.polygon(star, fill=(232, 72, 60, 255))
img.save(os.path.join(C, "Resources", "TomoSTAR.icns"), sizes=[(16, 16), (32, 32), (64, 64), (128, 128), (256, 256), (512, 512), (1024, 1024)])

# Zip with the Unix permissions in the entries, which Finder's Archive Utility and ditto restore.
z = os.path.join(out, f"TomoSTAR-{version}-macos.zip")
with zipfile.ZipFile(z, "w", zipfile.ZIP_DEFLATED) as zf:
    now = time.localtime()[:6]
    for root, dirs, files in os.walk(app):
        names = ([""] if root == app else []) + sorted(dirs) + sorted(files)
        for name in names:
            full = os.path.join(root, name) if name else root
            rel = os.path.relpath(full, out) + ("/" if os.path.isdir(full) else "")
            info = zipfile.ZipInfo(rel, date_time=now)
            info.create_system = 3
            mode = os.stat(full).st_mode
            info.external_attr = (mode & 0xFFFF) << 16 | (0x10 if os.path.isdir(full) else 0)
            if os.path.isdir(full):
                zf.writestr(info, "")
            else:
                info.compress_type = zipfile.ZIP_DEFLATED
                with open(full, "rb") as f:
                    zf.writestr(info, f.read())
print(f"{app}\n{z} ({os.path.getsize(z) / 1e6:.0f} MB)")
