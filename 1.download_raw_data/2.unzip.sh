#!/usr/bin/env bash
set -euo pipefail

#######################################
# CONFIGURATION
#######################################

DIRS=(
  "/data/srlab/AMP_collab/data/early_disease_synovium/xenium/Xenium_RA-SYN/Level_1_Xenium_Output"
  "/data/srlab/AMP_collab/data/early_disease_synovium/xenium/Xenium_PsD-SYN/Level_1_Xenium_Output"
  "/data/srlab/AMP_collab/data/early_disease_synovium/xenium/Xenium_CTRL-SYN/Level_1_Xenium_Output"
)


MASK="*.zip"

#######################################
# LOGIC
#######################################

for ROOT in "${DIRS[@]}"; do
  echo "Processing root directory: $ROOT"

  find "$ROOT" -type f -name "$MASK" | while read -r ZIPFILE; do
    echo "  Found: $ZIPFILE"
    echo "    Unzipping into: $ROOT"
    7zzs x -y "$ZIPFILE" -o"$ROOT"
    echo "========================="
  done

  echo "========================="
  echo "========================="
done

echo
read -r -p "Delete zip files in $ROOT? [y/N] " RESP
if [[ "$RESP" =~ ^[Yy]$ ]]; then
  for ROOT in "${DIRS[@]}"; do  
    echo "deleting files in $ROOT"
    rm -v "$ROOT"/*.zip
  done
else
  echo "Skipping deletion in $ROOT"
fi

echo "Warning: one of the zip files in the RA group also seems to contain extra xenium output that is"
echo "unnecessary and gets unzipped outside of that slide's folder. It is safe to delete these files."
