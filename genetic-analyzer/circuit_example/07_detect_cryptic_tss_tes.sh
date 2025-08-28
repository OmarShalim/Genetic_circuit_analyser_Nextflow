#!/usr/bin/env bash
# 07_detect_crytpic_tss_tes.sh
# Wrapper to call ../bin/detect_cryptic_tss_tes.py and write one TSV (known + cryptic)
#
# Usage:
#   ./07_detect_crytpic_tss_tes.sh [DATA_DIR] [GFF_FILE_OR_DIR] [WINDOWS] [MARGIN_BP]
#
# Defaults:
#   DATA_DIR  = current directory (scanned recursively: results/tube_1, tube_2, ...)
#   GFF_FILE  = auto-detect in DATA_DIR/data/gff (then DATA_DIR/gff, then DATA_DIR)
#   WINDOWS   = 3
#   MARGIN_BP = 1500
#
# Override Python with $PYTHON and script path with $SCRIPT_PATH if needed.

set -euo pipefail

DATA_DIR="${1:-.}"
GFF_INPUT="${2:-}"
WINDOWS="${3:-3}"
MARGIN="${4:-1500}"

PYTHON_BIN="${PYTHON:-python3}"

# Resolve this script's directory and default Python script location
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
SCRIPT_PATH_DEFAULT="${SCRIPT_DIR}/../bin/detect_cryptic_tss_tes.py"
SCRIPT_PATH="${SCRIPT_PATH:-$SCRIPT_PATH_DEFAULT}"

# Ensure DATA_DIR exists
if [[ ! -d "$DATA_DIR" ]]; then
  echo "ERROR: DATA_DIR not found: $DATA_DIR" >&2
  exit 1
fi

resolve_gff() {
  local datadir="$1"
  local provided="$2"

  if [[ -n "$provided" ]]; then
    if [[ -d "$provided" ]]; then
      local cand
      cand=$(ls -1 "$provided"/*.gff "$provided"/*.gff3 "$provided"/*.gtf 2>/dev/null | head -n 1 || true)
      if [[ -n "${cand:-}" ]]; then
        echo "$cand"; return 0
      else
        echo "ERROR: No .gff/.gff3/.gtf found in directory: $provided" >&2
        return 1
      fi
    elif [[ -f "$provided" ]]; then
      echo "$provided"; return 0
    else
      echo "ERROR: GFF path not found: $provided" >&2
      return 1
    fi
  fi

  # Prefer DATA_DIR/data/gff (your layout), then DATA_DIR/gff, then DATA_DIR
  local search_dirs=("$datadir/data/gff" "$datadir/gff" "$datadir")
  for d in "${search_dirs[@]}"; do
    if [[ -d "$d" ]]; then
      local cand
      cand=$(ls -1 "$d"/*.gff "$d"/*.gff3 "$d"/*.gtf 2>/dev/null | head -n 1 || true)
      if [[ -n "${cand:-}" ]]; then
        echo "$cand"; return 0
      fi
    fi
  done

  echo "ERROR: Could not find a GFF file. Provide one explicitly or place it under: $datadir/data/gff" >&2
  return 1
}

GFF_FILE="$(resolve_gff "$DATA_DIR" "$GFF_INPUT")"

# Ensure the Python script exists
if [[ ! -f "$SCRIPT_PATH" ]]; then
  if [[ -f "${DATA_DIR}/detect_cryptic_tss_tes.py" ]]; then
    SCRIPT_PATH="${DATA_DIR}/detect_cryptic_tss_tes.py"
  else
    echo "ERROR: detect_cryptic_tss_tes.py not found. Expected at $SCRIPT_PATH" >&2
    exit 1
  fi
fi

echo "DATA_DIR:  $DATA_DIR"
echo "GFF_FILE:  $GFF_FILE"
echo "WINDOWS:   $WINDOWS"
echo "MARGIN_BP: $MARGIN"
echo "PY_SCRIPT: $SCRIPT_PATH"
echo "Running: $PYTHON_BIN $SCRIPT_PATH --base-dir \"$DATA_DIR\" --gff \"$GFF_FILE\" --windows $WINDOWS --margin $MARGIN"

"$PYTHON_BIN" "$SCRIPT_PATH" --base-dir "$DATA_DIR" --gff "$GFF_FILE" --windows "$WINDOWS" --margin "$MARGIN"

