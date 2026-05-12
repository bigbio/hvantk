#!/usr/bin/env bash
# Verify Hail and Java are available for skill-driven sessions.

if ! command -v java >/dev/null 2>&1; then
  echo "[hvantk-skills] WARNING: java not found. Install JDK 11 via your OS package manager, then restart the shell."
  exit 0
fi

if ! python -c "import hail" >/dev/null 2>&1; then
  echo "[hvantk-skills] WARNING: hail not importable. Run 'poetry install'."
  exit 0
fi

echo "[hvantk-skills] hail + java available. Resource skills can run."
