#!/usr/bin/env bash
#
# Put "MOF SBU Viewer" in the applications menu, so starting it is a click and not a
# command.  Run once per machine:
#
#     ./launch/install-desktop-entry.sh
#     ./launch/install-desktop-entry.sh --gpu      # that shortcut uses the GPU
#
# Writes one file to ~/.local/share/applications and nothing else.  Remove it with
#     rm ~/.local/share/applications/mofsbu-viewer.desktop
set -euo pipefail

SOURCE=${BASH_SOURCE[0]}
while [ -L "$SOURCE" ]; do
  DIR=$(cd -P "$(dirname "$SOURCE")" && pwd); SOURCE=$(readlink "$SOURCE")
  [[ $SOURCE != /* ]] && SOURCE=$DIR/$SOURCE
done
HERE=$(cd -P "$(dirname "$SOURCE")" && pwd)
REPO=$(cd -P "$HERE/.." && pwd)

EXTRA=""
NAME="MOF SBU Viewer"
SLUG="mofsbu-viewer"
if [ "${1:-}" = "--gpu" ]; then
  EXTRA=" --gpu"; NAME="MOF SBU Viewer (GPU)"; SLUG="mofsbu-viewer-gpu"
fi

DEST="$HOME/.local/share/applications"
mkdir -p "$DEST"

# Terminal=true on purpose.  The server runs in the foreground for as long as you are
# using it, and the window is where its output goes — including the line that says which
# database it opened.  A silent background process that has to be hunted down in a
# process list to stop is not an improvement.
cat > "$DEST/$SLUG.desktop" <<EOF
[Desktop Entry]
Type=Application
Name=$NAME
Comment=Browse the mofsbu structure registry and queue builds
Exec=$HERE/mofsbu.sh$EXTRA
Path=$REPO
Terminal=true
Categories=Science;Chemistry;
EOF

chmod +x "$DEST/$SLUG.desktop"
command -v update-desktop-database >/dev/null 2>&1 && \
  update-desktop-database "$DEST" >/dev/null 2>&1 || true

echo "installed: $DEST/$SLUG.desktop"
echo "It should appear in the applications menu as \"$NAME\"."
echo "You can also drag that file to the desktop to get a double-clickable icon."
