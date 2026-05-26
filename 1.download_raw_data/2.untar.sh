#!/usr/bin/env bash
set -euo pipefail

#######################################
# CONFIGURATION
#######################################

DIRS=(
   "/data/srlab/AMP_collab/data/early_disease_synovium/xenium/Xenium_RA-SYN-ARBITRATE/Level_1_Xenium_Output/pre_treatment/"
   "/data/srlab/AMP_collab/data/early_disease_synovium/xenium/Xenium_RA-SYN-ARBITRATE/Level_1_Xenium_Output/post_treatment/"
)

MASK="*.tar.gz"

#######################################
# LOGIC
#######################################

for ROOT in "${DIRS[@]}"; do
  echo "Processing root directory: $ROOT"

  find "$ROOT" -type f -name "$MASK" | while read -r ZIPFILE; do
    echo "  Found: $ZIPFILE"
    echo "    Unzipping into: $ROOT"
    pigz -kdc $ZIPFILE | tar xf - -C $ROOT
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
