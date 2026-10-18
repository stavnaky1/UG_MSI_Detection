#!/usr/bin/env bash
# Install check: run the synthetic toy example and compare to expected_output.csv.
# toy_msi carries INFO/GNOMAD_AF; toy_mss does not, so this also exercises --gnomad_tsv.
set -euo pipefail
cd "$(dirname "$0")"
PY="${PYTHON:-python}"

rm -rf toy_out
"$PY" hmer_indel_summary.py \
    --bed_for_intersect toy_regions.bed \
    --path_with_vcfs    vcfs \
    --output_dir        toy_out \
    --gnomad_tsv        toy_gnomad_af.tsv.gz

"$PY" - <<'PYEOF'
import pandas as pd
obs = pd.read_csv("toy_out/vcfs_toy_regions_hmer_indel_summary.csv")
exp = pd.read_csv("expected_output.csv")
pd.testing.assert_frame_equal(obs, exp, check_exact=False, rtol=1e-9)
print("\nOK: output matches expected_output.csv")
PYEOF
