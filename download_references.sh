#!/usr/bin/env bash
# Download the reference BED files from Zenodo (doi:10.5281/zenodo.23160540)
# and verify their md5 checksums.
#
# Usage: bash scripts/download_references.sh [output_dir]   (default: refs/)
set -euo pipefail

OUT="${1:-refs}"
RECORD="https://zenodo.org/records/23160540/files"
mkdir -p "$OUT"

md5_of() {
    if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | cut -d' ' -f1
    else md5 -q "$1"; fi   # macOS
}

fetch() {  # $1 = file name, $2 = expected md5
    local f="$1" expected="$2" dest="$OUT/$1"
    if [[ -f "$dest" && "$(md5_of "$dest")" == "$expected" ]]; then
        echo "OK (cached)  $f"; return
    fi
    echo "Downloading  $f"
    curl -fL --retry 3 -o "$dest.part" "$RECORD/$f?download=1"
    mv "$dest.part" "$dest"
    local got; got="$(md5_of "$dest")"
    if [[ "$got" != "$expected" ]]; then
        echo "ERROR: md5 mismatch for $f (expected $expected, got $got)" >&2
        exit 1
    fi
    echo "OK           $f"
}

fetch ug_hcr.v1.4.bed \
      760b72fc727a475a27e47d3c396e0d7a
fetch hg38_exome_v2.0.2_targets_sorted_validated.annotated_0.padded_100.bed \
      0c0f17a6732c82bfeda9795fe839e5c4

echo "Reference BEDs ready in $OUT/"
