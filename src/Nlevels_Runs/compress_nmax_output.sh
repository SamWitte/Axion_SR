#!/bin/bash
# Compress one completed Nmax run (all five .dat files) to NL_output_Nmax_<N>.tar.gz
# Usage: compress_nmax_output.sh <outdir> <Nmax>
# Removes uncompressed .dat only after the archive passes integrity check.
set -euo pipefail

OUTDIR="${1:?outdir required}"
NMAX="${2:?Nmax required}"
ARCHIVE="NL_output_Nmax_${NMAX}.tar.gz"

cd "$OUTDIR"

if [[ -f "$ARCHIVE" && "${FORCE:-0}" != "1" ]]; then
    echo "Already compressed: ${OUTDIR}/${ARCHIVE}"
    exit 0
fi

prefixes=(Time Spin States Modes MassBH)
files=()
for p in "${prefixes[@]}"; do
    matches=( "${p}"_*Nmax_"${NMAX}".dat )
    if [[ ! -f "${matches[0]}" ]]; then
        echo "Skip compress (incomplete): missing ${p} for Nmax=${NMAX} in ${OUTDIR}"
        exit 0
    fi
    # Empty States/Modes are valid after pruning removed all inactive modes.
    files+=( "${matches[0]}" )
done

echo "Compressing Nmax=${NMAX} in ${OUTDIR} (${#files[@]} files) ..."
if command -v pigz >/dev/null 2>&1; then
    tar -cf - "${files[@]}" | pigz -9 > "${ARCHIVE}.tmp"
    pigz -t "${ARCHIVE}.tmp"
else
    tar -czf "${ARCHIVE}.tmp" "${files[@]}"
    gzip -t "${ARCHIVE}.tmp"
fi
mv "${ARCHIVE}.tmp" "$ARCHIVE"
rm -f "${files[@]}"
echo "Done: ${OUTDIR}/${ARCHIVE}"
