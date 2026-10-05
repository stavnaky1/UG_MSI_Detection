#!/usr/bin/env python3
"""
Created on Mon Jan 12 10:42:46 2026

@author: stavnaky

Usage:
  python hmer_indel_summary.py \
      --bed_for_intersect /storage/bfe_maruvka/stavnbfe/UG_MSI/hg38_exome_v2.0.2_targets_sorted_validated.annotated_0.padded_100.bed \
      --path_with_vcfs /storage/bfe_maruvka/stavnbfe/UG_MSI/dilution_vcfs \
      --output_dir /storage/bfe_maruvka/stavnbfe/UG_MSI/dilution_vcfs/_filtered_cache_exome_17062026
"""

# ======================================================
# Imports
# ======================================================
import os
import glob
import shlex
import argparse
import subprocess
from io import StringIO

import numpy as np
import pandas as pd


# Output-name tag per BED stem, so outputs keep the old
# <vcf dir name>_<wgs|wxs>_hmer_indel_summary.csv names.
# Unknown BEDs fall back to the BED stem.
BED_TAGS = {
    "ug_hcr.v1.4": "wgs",
    "hg38_exome_v2.0.2_targets_sorted_validated.annotated_0.padded_100": "wxs",
}


# ======================================================
# Command-line arguments
# ======================================================
def parse_args():
    p = argparse.ArgumentParser(
        description="Summarize homopolymer indels in somatic VCFs within a BED territory.")
    p.add_argument("-b", "--bed_for_intersect", required=True,
                   help="BED to intersect with (whole genome or 100bp-padded exome).")
    p.add_argument("-v", "--path_with_vcfs", required=True,
                   help="Directory containing the input *.vcf.gz files.")
    p.add_argument("-o", "--output_dir", required=True,
                   help="Output directory for prefiltered VCFs and the summary CSV "
                        "(created if missing; must differ from --path_with_vcfs).")
    return p.parse_args()


# ======================================================
# MAIN FUNCTION (top)
# ======================================================
def main():
    global output_dir, limit_to_at, min_qual, max_indel_length
    global bed_stem, bed_for_intersect_global, bcftools_threads
    global gnomad_threshold, af_thresh, min_indel_length

    args = parse_args()

    # --------------------------------------------------
    # Inputs (command line)
    # --------------------------------------------------
    bed_for_intersect = args.bed_for_intersect
    path_with_vcfs = os.path.normpath(args.path_with_vcfs)
    output_dir = os.path.normpath(args.output_dir)
    assert os.path.realpath(output_dir) != os.path.realpath(path_with_vcfs), \
        "output_dir must differ from path_with_vcfs (filtered VCFs would be re-read as inputs)"

    bed_for_intersect_global = bed_for_intersect
    bed_stem = os.path.splitext(os.path.basename(bed_for_intersect))[0]

    out_tag = BED_TAGS.get(bed_stem, bed_stem)
    out_path = os.path.join(
        output_dir,
        f"{os.path.basename(path_with_vcfs)}_{out_tag}_hmer_indel_summary.csv")

    # --------------------------------------------------
    # Hard-coded parameters
    # --------------------------------------------------
    limit_to_at = False
    af_thresh = 0.5
    min_qual = 16

    gnomad_threshold = 0.01
    min_indel_length = 10
    max_indel_length = 15

    bcftools_threads = 8

    print("bed_for_intersect:     ", bed_for_intersect)
    print("path_with_vcfs:        ", path_with_vcfs)
    print("output_dir:            ", output_dir)
    print("out_path:              ", out_path)

    # --------------------------------------------------
    # VCF discovery
    # --------------------------------------------------
    # vcf_paths = sorted(glob.glob(os.path.join(path_with_vcfs, '_filtered_cache/*.vcf.gz')))
    vcf_paths = sorted(glob.glob(os.path.join(path_with_vcfs, '*.vcf.gz')))
    vcf_paths = [p for p in vcf_paths if not p.endswith('.tbi') and not p.endswith('.csi')]

    print('VCFs found:', len(vcf_paths))
    print('Example VCFs:', [os.path.basename(p) for p in vcf_paths[:5]])
    assert os.path.exists(bed_for_intersect), f"bed_for_intersect not found: {bed_for_intersect}"

    # --------------------------------------------------
    # Step 1: Prefilter/cached VCFs
    # --------------------------------------------------
    os.makedirs(output_dir, exist_ok=True)

    filtered_vcf_paths = []
    for i, vcf in enumerate(vcf_paths, start=1):
        if i % 10 == 1 or i == len(vcf_paths):
            print(f"Prefiltering {i}/{len(vcf_paths)}: {os.path.basename(vcf)}")
        filtered_vcf_paths.append(_ensure_filtered_vcf(vcf))
    # filtered_vcf_paths = vcf_paths

    print("Filtered VCFs ready:", len(filtered_vcf_paths))

    # --------------------------------------------------
    # Step 2: Measurements
    # --------------------------------------------------
    rows = []

    for i, vcf_path in enumerate(filtered_vcf_paths, start=1):
        vcf_name = os.path.basename(vcf_path)

        # if i % 10 == 1 or i == len(filtered_vcf_paths):
        print(f"Measuring {i}/{len(filtered_vcf_paths)}: {vcf_name}")

        try:
            raw_df = _bcftools_query_df(vcf_path)
        except:
            continue
        n_in_bed = len(raw_df)
        df = raw_df.copy()

        # Apply filters
        df = df[df["QUAL"].fillna(-np.inf) >= min_qual]
        df = df[~df["ALT"].astype(str).str.contains(",", regex=False)]
        df = df[df["X_HIL"] <= max_indel_length]

        if limit_to_at:
            df = df[df["X_HIN"].isin(["A", "T"])]
        df["GNOMAD_AF"] = pd.to_numeric(df["GNOMAD_AF"], errors="coerce")
        df = df[(df["ID"] == ".") & ((df["GNOMAD_AF"] < gnomad_threshold) | (df["GNOMAD_AF"].isna()))]
        df = df[df["VAF"].fillna(np.inf) <= af_thresh]

        n_after_filters = len(df)

        # Empty case
        if n_after_filters == 0:
            rows.append(_empty_row(vcf_name, vcf_path, n_in_bed))
            continue
        row = _measure_row(df, vcf_name, vcf_path, n_in_bed, n_after_filters)
        rows.append(row)
        print(row)
    # --------------------------------------------------
    # Final DataFrame
    # --------------------------------------------------
    result_df = pd.DataFrame(rows)
    # print(result_df.head())
    result_df.to_csv(out_path, index=False)
    print(f"Saved summary → {out_path}")

# ======================================================
# HELPER FUNCTIONS (bottom)
# ======================================================

def _ensure_gnomad(in_vcf):
    """
    Ensure that the VCF has GNOMAD_AF field.
    If missing, annotate using bcftools annotate.
    Returns path to VCF with GNOMAD_AF guaranteed.
    """

    # output path (same folder as input)
    annotated_vcf = in_vcf.replace(".vcf.gz", ".gnomad_annotated.vcf.gz")
    annotated_tbi = annotated_vcf + ".tbi"

    # 1. Check if GNOMAD_AF exists in header
    check_cmd = f"bcftools view -h {shlex.quote(in_vcf)} | grep -q 'ID=GNOMAD_AF'"
    has_af = subprocess.call(["bash", "-lc", check_cmd]) == 0

    if has_af:
        return in_vcf

    if os.path.exists(annotated_vcf) and os.path.exists(annotated_tbi):
        return annotated_vcf

    print(f"Annotating GNOMAD_AF → {os.path.basename(annotated_vcf)}")

    gnomad_tsv = "/storage/bfe_maruvka/References/GRCh38_chromosomes_vcf/hg38_snps_all_chromosomes_with_af.tsv.gz"

    annotate_cmd = (
        "set -euo pipefail; "
        f"bcftools annotate "
        f"-a {shlex.quote(gnomad_tsv)} "
        "-c CHROM,POS,REF,ALT,INFO/GNOMAD_AF "
        "-h <(echo '##INFO=<ID=GNOMAD_AF,Number=A,Type=Float,Description=\"Gnomad Allele Frequency\">') "
        f"-O z -o {shlex.quote(annotated_vcf)} "
        f"{shlex.quote(in_vcf)}; "
        f"bcftools index -t {shlex.quote(annotated_vcf)}"
    )

    subprocess.check_call(["bash", "-lc", annotate_cmd])

    return annotated_vcf
def _empty_row(vcf_name, vcf_path, n_in_bed):
    return {
        "vcf": vcf_name,
        "vcf_path": vcf_path,
        "n_records_in_bed": n_in_bed,
        "n_records_after_filters": 0,
        "Total indels": 0,
        "h-indels": 0,
        "% h-indels": np.nan,
        "1nt del": 0,
        "% 1nt del": np.nan,
        "2nt del": 0,
        "1nt ins": 0,
        "Long h-indels": 0,
        "% Long h-indels in all indels": np.nan,
        "1nt del in long hmers": 0,
        "% 1nt del in long hmers": np.nan,
        "% 1nt del in long hmers of all indels": np.nan,
        "1 nt del / 1 nt ins": np.nan,
        "2 nt del / 1 nt ins": np.nan,
    }


def _measure_row(df, vcf_name, vcf_path, n_in_bed, n_after_filters):

    vt = df["VARIANT_TYPE"].astype(str)

    n_hindel = int((vt == "h-indel").sum())
    n_nonh_indel = int((vt == "non-h-indel").sum())
    total_indels = n_hindel + n_nonh_indel

    idx_hindels = (vt == "h-indel").to_numpy(bool)

    ref = df["REF"].astype(str).to_numpy()
    alt = df["ALT"].astype(str).to_numpy()

    one_nt_del = np.asarray([(len(r) - len(a)) == 1 for r, a in zip(ref, alt)], dtype=bool)
    two_nt_del = np.asarray([(len(r) - len(a)) == 2 for r, a in zip(ref, alt)], dtype=bool)
    one_nt_ins = np.asarray([(len(a) - len(r)) == 1 for r, a in zip(ref, alt)], dtype=bool)

    x_hil = df["X_HIL"].to_numpy()
    idx_long_hindels = ((x_hil >= min_indel_length) & (x_hil <= max_indel_length) & idx_hindels).astype(bool)

    n_1nt_del = int((one_nt_del & idx_hindels).sum())
    n_2nt_del = int((two_nt_del & idx_hindels).sum())
    n_1nt_ins = int((one_nt_ins & idx_hindels).sum())
    n_long_hindels = int(idx_long_hindels.sum())
    n_1nt_del_in_long = int((one_nt_del & idx_long_hindels).sum())

    return {
        "vcf": vcf_name,
        "vcf_path": vcf_path,
        "n_records_in_bed": n_in_bed,
        "n_records_after_filters": n_after_filters,
        "Total indels": total_indels,
        "h-indels": n_hindel,
        "% h-indels": _safe_pct(n_hindel, total_indels),
        "1nt del": n_1nt_del,
        "% 1nt del": _safe_pct(n_1nt_del, total_indels),
        "2nt del": n_2nt_del,
        "1nt ins": n_1nt_ins,
        "Long h-indels": n_long_hindels,
        "% Long h-indels in all indels": _safe_pct(n_long_hindels, total_indels),
        "1nt del in long hmers": n_1nt_del_in_long,
        "% 1nt del in long hmers": _safe_pct(n_1nt_del_in_long, n_long_hindels),
        "% 1nt del in long hmers of all indels": _safe_pct(n_1nt_del_in_long, total_indels),
        "1 nt del / 1 nt ins": (n_1nt_del / n_1nt_ins) if n_1nt_ins > 0 else np.nan,
        "2 nt del / 1 nt ins": (n_2nt_del / n_1nt_ins) if n_1nt_ins > 0 else np.nan,
    }


def _safe_pct(numer, denom):
    if denom in (None, 0) or (isinstance(denom, float) and np.isnan(denom)):
        return np.nan
    return numer / denom * 100


def _filtered_vcf_path(in_vcf):
    base = os.path.basename(in_vcf)
    stem = base[:-7] if base.endswith('.vcf.gz') else base
    at_tag = 'AT' if limit_to_at else 'all'
    return os.path.join(
        output_dir,
        f"{stem}.flt_{bed_stem}_PASS_indels_q{min_qual}_x{max_indel_length}_{at_tag}.vcf.gz",
    )


def _ensure_filtered_vcf(in_vcf):
    out_vcf = _filtered_vcf_path(in_vcf)
    out_tbi = out_vcf + ".tbi"

    if os.path.exists(out_vcf):
        if not os.path.exists(out_tbi):
            subprocess.check_call([
                "bash", "-lc",
                f"bcftools index -t --threads {bcftools_threads} {shlex.quote(out_vcf)}"
            ])
        return _ensure_gnomad(out_vcf)

    expr_parts = [
    f"QUAL>={min_qual}",
    "ID=='.'",
    f"INFO/X_HIL<={max_indel_length}",
    'ALT!~","'
    ]

    if limit_to_at:
        expr_parts.append("(INFO/X_HIN=='A' || INFO/X_HIN=='T')")

    expr = " && ".join(expr_parts)

    cmd = (
        "set -euo pipefail; "
        f"bcftools view --threads {bcftools_threads} -R {shlex.quote(bed_for_intersect_global)} "
        f"-f PASS -v indels -i {shlex.quote(expr)} -Oz -o {shlex.quote(out_vcf)} {shlex.quote(in_vcf)}; "
        f"bcftools index -t --threads {bcftools_threads} {shlex.quote(out_vcf)}"
    )

    print(f"Filtering -> {os.path.basename(out_vcf)}")
    subprocess.check_call(["bash", "-lc", cmd])

    return _ensure_gnomad(out_vcf)


def _bcftools_query_df(vcf_path):
    fmt = (
        "%CHROM\t%POS\t%REF\t%ALT\t%ID\t%QUAL\t"
        "%INFO/VARIANT_TYPE\t%INFO/X_HIL\t%INFO/X_HIN\t%INFO/X_IC\t%INFO/GNOMAD_AF\t[%VAF]\n"
    )

    cmd = f"bcftools query -f {shlex.quote(fmt)} {shlex.quote(vcf_path)}"
    out = subprocess.check_output(["bash", "-lc", cmd], text=True)

    cols = [
        "CHROM", "POS", "REF", "ALT", "ID", "QUAL",
        "VARIANT_TYPE", "X_HIL", "X_HIN", "X_IC",
        "GNOMAD_AF", "VAF",
    ]

    if not out.strip():
        return pd.DataFrame(columns=cols)

    df = pd.read_csv(StringIO(out), sep="\t", names=cols)

    df["POS"] = pd.to_numeric(df["POS"], errors="coerce").astype("Int64")
    df["QUAL"] = pd.to_numeric(df["QUAL"], errors="coerce")
    df["VAF"] = pd.to_numeric(df["VAF"], errors="coerce")

    df["X_HIL"] = pd.to_numeric(df["X_HIL"], errors="coerce").fillna(0).astype(int)

    df["GNOMAD_AF"] = df["GNOMAD_AF"].replace({"-": np.nan, ".": np.nan, "": np.nan})
    df["GNOMAD_AF"] = pd.to_numeric(df["GNOMAD_AF"], errors="coerce")

    for c in ["ID", "VARIANT_TYPE", "X_HIN", "X_IC"]:
        df[c] = df[c].fillna(".")

    return df


# ======================================================
# Run main if executed directly
# ======================================================
if __name__ == "__main__":
    main()