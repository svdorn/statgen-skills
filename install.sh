#!/usr/bin/env bash
# install.sh - symlink each skill in this repo into ~/.claude/skills/
#
# Usage:
#   ./install.sh          # symlink all skills
#   ./install.sh <name>   # symlink a single named skill
#
# Idempotent: skips existing symlinks pointing at the right target.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CLAUDE_SKILLS_DIR="${CLAUDE_SKILLS_DIR:-$HOME/.claude/skills}"
mkdir -p "$CLAUDE_SKILLS_DIR"

link_skill() {
  local name="$1"
  local src="$REPO_ROOT/$name"
  local dst="$CLAUDE_SKILLS_DIR/$name"

  if [[ ! -f "$src/SKILL.md" ]]; then
    echo "SKIP: $name (no SKILL.md)"; return 0
  fi
  if [[ -L "$dst" && "$(readlink "$dst")" == "$src" ]]; then
    echo "OK:   $name (already symlinked)"; return 0
  fi
  if [[ -e "$dst" ]]; then
    echo "WARN: $dst exists and is not the expected symlink; skipping" >&2
    return 1
  fi
  ln -s "$src" "$dst"
  echo "LINK: $name -> $src"
}

if [[ $# -gt 0 ]]; then
  link_skill "$1"
else
  for d in "$REPO_ROOT"/*/; do
    name="$(basename "$d")"
    [[ "$name" == "scripts" || "$name" == "references" || "$name" == "assets" ]] && continue
    link_skill "$name"
  done
fi

echo "Done. Restart your agent client to pick up new skills."
