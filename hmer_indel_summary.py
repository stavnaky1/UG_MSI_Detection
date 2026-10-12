#!/usr/bin/env python3
"""
hmer_indel_summary.py

Homopolymer-indel (h-indel) summary statistics for microsatellite instability
(MSI) analysis of Ultima Genomics somatic VCFs.

For every *.vcf.gz in an input directory, the script:
  1. prefilters the VCF with bcftools to a BED territory (whole-genome HCR or
     padded exome), PASS indels, QUAL >= --min_qual, no known ID, homopolymer
     length <= --max_hmer_length (cached in --output_dir);
  2. adds INFO/GNOMAD_AF from --gnomad_tsv if the VCF does not already carry it;
  3. removes multiallelic records, gnomAD-common variants and variants with
     VAF > --max_vaf (or missing VAF);
  4. counts h-indels, 1-nt deletions, long-homopolymer 1-nt deletions, etc.,
     and writes one row per VCF to a summary CSV.

Example:
  python hmer_indel_summary.py \\
      --bed_for_intersect refs/ug_hcr.v1.4.bed \\
      --path_with_vcfs    my_vcfs/ \\
      --output_dir        my_vcfs_out/ \\
      --gnomad_tsv        refs/gnomad_af.tsv.gz

See README.md for input requirements and output column definitions.
"""

import argparse
import glob
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from io import StringIO
from typing import Optional

import numpy as np
import pandas as pd

__version__ = "1.0.0"

# Output-name tag per BED stem: <vcf dir name>_<tag>_hmer_indel_summary.csv.
# BEDs not listed here fall back to their file stem.
BED_TAGS = {
    "ug_hcr.v1.4": "wgs",
    "hg38_exome_v2.0.2_targets_sorted_validated.annotated_0.padded_100": "wxs",
}

QUERY_FORMAT = (
    "%CHROM\t%POS\t%REF\t%ALT\t%ID\t%QUAL\t"
    "%INFO/VARIANT_TYPE\t%INFO/X_HIL\t%INFO/X_HIN\t%INFO/GNOMAD_AF\t[%VAF]\n"
)
QUERY_COLS = [
    "CHROM", "POS", "REF", "ALT", "ID", "QUAL",
    "VARIANT_TYPE", "X_HIL", "X_HIN", "GNOMAD_AF", "VAF",
]
REQUIRED_TAGS = [("INFO", "VARIANT_TYPE"), ("INFO", "X_HIL"), ("INFO", "X_HIN"),
                 ("FORMAT", "VAF")]
GNOMAD_HEADER_LINE = (
    '##INFO=<ID=GNOMAD_AF,Number=A,Type=Float,Description="Gnomad Allele Frequency">\n'
)


# ======================================================
# Parameters
# ======================================================
@dataclass(frozen=True)
class Params:
    bed: str
    output_dir: str
    gnomad_tsv: Optional[str]
    min_qual: float
    max_hmer_length: int
    min_long_hmer_length: int
    gnomad_af_threshold: float
    max_vaf: float
    at_only: bool
    threads: int

    @property
    def bed_stem(self):
        return os.path.splitext(os.path.basename(self.bed))[0]


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Summarize homopolymer indels in Ultima Genomics somatic VCFs "
                    "within a BED territory. Defaults reproduce the published analysis.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    req = p.add_argument_group("required")
    req.add_argument("-b", "--bed_for_intersect", required=True,
                     help="BED territory (e.g. ug_hcr.v1.4.bed for WGS, or the "
                          "100bp-padded exome BED for WXS).")
    req.add_argument("-v", "--path_with_vcfs", required=True,
                     help="Directory containing the input *.vcf.gz files "
                          "(bgzipped, indexed, single-sample).")
    req.add_argument("-o", "--output_dir", required=True,
                     help="Output directory for cached filtered VCFs and the summary "
                          "CSV (created if missing; must differ from --path_with_vcfs).")

    p.add_argument("-g", "--gnomad_tsv", default=None,
                   help="bgzipped, tabix-indexed TSV (CHROM POS REF ALT AF) used to add "
                        "INFO/GNOMAD_AF to VCFs that lack it. Required only for such VCFs.")

    flt = p.add_argument_group("filters")
    flt.add_argument("--min_qual", type=float, default=16,
                     help="Minimum QUAL.")
    flt.add_argument("--max_hmer_length", type=int, default=15,
                     help="Maximum homopolymer length (INFO/X_HIL); longer records "
                          "are removed.")
    flt.add_argument("--min_long_hmer_length", type=int, default=10,
                     help="Minimum INFO/X_HIL for an h-indel to count as 'long'.")
    flt.add_argument("--gnomad_af_threshold", type=float, default=0.01,
                     help="Keep variants with GNOMAD_AF below this value (or absent).")
    flt.add_argument("--max_vaf", type=float, default=0.5,
                     help="Keep variants with FORMAT/VAF <= this value "
                          "(variants with missing VAF are removed).")
    flt.add_argument("--at_only", action="store_true",
                     help="Restrict to A/T homopolymers (INFO/X_HIN in {A,T}).")

    p.add_argument("-t", "--threads", type=int, default=8,
                   help="bcftools compression/decompression threads.")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    args = p.parse_args(argv)

    if os.path.realpath(args.output_dir) == os.path.realpath(args.path_with_vcfs):
        p.error("--output_dir must differ from --path_with_vcfs "
                "(filtered VCFs would be re-read as inputs).")
    if not os.path.isfile(args.bed_for_intersect):
        p.error(f"BED not found: {args.bed_for_intersect}")
    if not os.path.isdir(args.path_with_vcfs):
        p.error(f"VCF directory not found: {args.path_with_vcfs}")
    if args.gnomad_tsv is not None and not os.path.isfile(args.gnomad_tsv):
        p.error(f"gnomAD TSV not found: {args.gnomad_tsv}")
    if args.min_long_hmer_length > args.max_hmer_length:
        p.error("--min_long_hmer_length cannot exceed --max_hmer_length.")
    return args


# ======================================================
# Main
# ======================================================
def main(argv=None):
    args = parse_args(argv)
    if shutil.which("bcftools") is None:
        sys.exit("ERROR: bcftools not found on PATH (see environment.yml).")

    path_with_vcfs = os.path.normpath(args.path_with_vcfs)
    params = Params(
        bed=args.bed_for_intersect,
        output_dir=os.path.normpath(args.output_dir),
        gnomad_tsv=args.gnomad_tsv,
        min_qual=args.min_qual,
        max_hmer_length=args.max_hmer_length,
        min_long_hmer_length=args.min_long_hmer_length,
        gnomad_af_threshold=args.gnomad_af_threshold,
        max_vaf=args.max_vaf,
        at_only=args.at_only,
        threads=args.threads,
    )

    out_tag = BED_TAGS.get(params.bed_stem, params.bed_stem)
    out_path = os.path.join(
        params.output_dir,
        f"{os.path.basename(path_with_vcfs)}_{out_tag}_hmer_indel_summary.csv")

    print(f"hmer_indel_summary.py v{__version__}")
    for k, v in [("bed", params.bed), ("path_with_vcfs", path_with_vcfs),
                 ("output_dir", params.output_dir), ("out_path", out_path),
                 ("gnomad_tsv", params.gnomad_tsv), ("min_qual", f"{params.min_qual:g}"),
                 ("max_hmer_length", params.max_hmer_length),
                 ("min_long_hmer_length", params.min_long_hmer_length),
                 ("gnomad_af_threshold", params.gnomad_af_threshold),
                 ("max_vaf", params.max_vaf), ("at_only", params.at_only)]:
        print(f"  {k:<21} {v}")

    vcf_paths = sorted(glob.glob(os.path.join(path_with_vcfs, "*.vcf.gz")))
    if not vcf_paths:
        sys.exit(f"ERROR: no *.vcf.gz files in {path_with_vcfs}")
    print(f"VCFs found: {len(vcf_paths)}")

    for vcf in vcf_paths:
        validate_input_vcf(vcf)
    if params.gnomad_tsv is None:
        lacking = [os.path.basename(v) for v in vcf_paths if not _has_info_tag(v, "GNOMAD_AF")]
        if lacking:
            sys.exit(f"ERROR: {len(lacking)} VCF(s) lack INFO/GNOMAD_AF (e.g. {lacking[0]}); "
                     "provide --gnomad_tsv.")
    os.makedirs(params.output_dir, exist_ok=True)

    # Step 1-2: prefilter (cached) and ensure GNOMAD_AF
    filtered_vcf_paths = []
    for i, vcf in enumerate(vcf_paths, start=1):
        if i % 10 == 1 or i == len(vcf_paths):
            print(f"Prefiltering {i}/{len(vcf_paths)}: {os.path.basename(vcf)}")
        filtered_vcf_paths.append(ensure_filtered_vcf(vcf, params))

    # Step 3-4: filter and measure
    rows, skipped = [], []
    for i, vcf_path in enumerate(filtered_vcf_paths, start=1):
        vcf_name = os.path.basename(vcf_path)
        print(f"Measuring {i}/{len(filtered_vcf_paths)}: {vcf_name}")
        try:
            raw_df = bcftools_query_df(vcf_path)
        except (subprocess.CalledProcessError, ValueError) as e:
            print(f"WARNING: skipping {vcf_name}: {e}", file=sys.stderr)
            skipped.append(vcf_name)
            continue
        df = apply_filters(raw_df, params)
        rows.append(measure_row(df, vcf_name, vcf_path, len(raw_df), params))

    columns = list(measure_row(pd.DataFrame(columns=QUERY_COLS), "", "", 0, params))
    pd.DataFrame(rows, columns=columns).to_csv(out_path, index=False)
    print(f"Saved summary ({len(rows)} rows) -> {out_path}")
    if skipped:
        print(f"WARNING: {len(skipped)} VCF(s) skipped and absent from the summary:",
              file=sys.stderr)
        for name in skipped:
            print(f"  {name}", file=sys.stderr)


# ======================================================
# bcftools helpers
# ======================================================
def _run(cmd):
    subprocess.run(cmd, check=True)


def _vcf_stem(path):
    base = os.path.basename(path)
    return base[:-len(".vcf.gz")] if base.endswith(".vcf.gz") else base


def filtered_vcf_path(in_vcf, params):
    at_tag = "AT" if params.at_only else "all"
    return os.path.join(
        params.output_dir,
        f"{_vcf_stem(in_vcf)}.flt_{params.bed_stem}_PASS_indels"
        f"_q{params.min_qual:g}_x{params.max_hmer_length}_{at_tag}.vcf.gz",
    )


def ensure_filtered_vcf(in_vcf, params):
    """Region/PASS/indel/QUAL/ID/hmer-length prefilter, cached by parameters."""
    out_vcf = filtered_vcf_path(in_vcf, params)

    if not os.path.exists(out_vcf):
        expr_parts = [
            f"QUAL>={params.min_qual:g}",
            "ID=='.'",
            f"INFO/X_HIL<={params.max_hmer_length}",
            'ALT!~","',
        ]
        if params.at_only:
            expr_parts.append("(INFO/X_HIN=='A' || INFO/X_HIN=='T')")
        expr = " && ".join(expr_parts)

        print(f"Filtering -> {os.path.basename(out_vcf)}")
        tmp = out_vcf + ".tmp"
        _run(["bcftools", "view", "--threads", str(params.threads),
              "-R", params.bed, "-f", "PASS", "-v", "indels", "-i", expr,
              "-Oz", "-o", tmp, in_vcf])
        os.replace(tmp, out_vcf)  # atomic: an interrupted run never leaves a partial cache

    if not os.path.exists(out_vcf + ".tbi"):
        _run(["bcftools", "index", "-t", "--threads", str(params.threads), out_vcf])

    return ensure_gnomad(out_vcf, params)


def _vcf_header(vcf):
    return subprocess.run(["bcftools", "view", "-h", vcf], check=True,
                          capture_output=True, text=True).stdout


def _has_tag(header, field, tag):
    return re.search(rf"^##{field}=<ID={re.escape(tag)},", header, re.MULTILINE) is not None


def _has_info_tag(vcf, tag):
    return _has_tag(_vcf_header(vcf), "INFO", tag)


def validate_input_vcf(vcf):
    """Fail early, naming the file, if a VCF cannot be processed."""
    name = os.path.basename(vcf)
    if not (os.path.exists(vcf + ".tbi") or os.path.exists(vcf + ".csi")):
        sys.exit(f"ERROR: {name} is not indexed (run: tabix -p vcf {name}).")
    header = _vcf_header(vcf)
    missing = [f"{field}/{tag}" for field, tag in REQUIRED_TAGS
               if not _has_tag(header, field, tag)]
    if missing:
        sys.exit(f"ERROR: {name} lacks required header fields: {', '.join(missing)}. "
                 "Inputs must be Ultima Genomics somatic VCFs carrying these annotations "
                 "(see README).")


def ensure_gnomad(in_vcf, params):
    """Return a path to a VCF carrying INFO/GNOMAD_AF, annotating if needed."""
    if _has_info_tag(in_vcf, "GNOMAD_AF"):
        return in_vcf

    annotated_vcf = in_vcf[:-len(".vcf.gz")] + ".gnomad_annotated.vcf.gz"

    if not os.path.exists(annotated_vcf):
        if params.gnomad_tsv is None:
            sys.exit(f"ERROR: {os.path.basename(in_vcf)} has no INFO/GNOMAD_AF; "
                     "provide --gnomad_tsv.")
        print(f"Annotating GNOMAD_AF -> {os.path.basename(annotated_vcf)}")
        tmp = annotated_vcf + ".tmp"
        with tempfile.NamedTemporaryFile("w", suffix=".hdr", delete=False) as fh:
            fh.write(GNOMAD_HEADER_LINE)
            hdr = fh.name
        try:
            _run(["bcftools", "annotate", "-a", params.gnomad_tsv,
                  "-c", "CHROM,POS,REF,ALT,INFO/GNOMAD_AF", "-h", hdr,
                  "-Oz", "-o", tmp, in_vcf])
        finally:
            os.remove(hdr)
        os.replace(tmp, annotated_vcf)

    if not os.path.exists(annotated_vcf + ".tbi"):
        _run(["bcftools", "index", "-t", annotated_vcf])

    return annotated_vcf


def bcftools_query_df(vcf_path):
    samples = subprocess.run(["bcftools", "query", "-l", vcf_path], check=True,
                             capture_output=True, text=True).stdout.split()
    if len(samples) != 1:
        raise ValueError(f"expected a single-sample VCF, found {len(samples)} samples")

    out = subprocess.run(["bcftools", "query", "-f", QUERY_FORMAT, vcf_path],
                         check=True, capture_output=True, text=True).stdout
    if not out.strip():
        return pd.DataFrame(columns=QUERY_COLS)

    df = pd.read_csv(StringIO(out), sep="\t", names=QUERY_COLS)

    df["POS"] = pd.to_numeric(df["POS"], errors="coerce").astype("Int64")
    df["QUAL"] = pd.to_numeric(df["QUAL"], errors="coerce")
    df["VAF"] = pd.to_numeric(df["VAF"], errors="coerce")
    df["X_HIL"] = pd.to_numeric(df["X_HIL"], errors="coerce").fillna(0).astype(int)
    df["GNOMAD_AF"] = pd.to_numeric(
        df["GNOMAD_AF"].replace({"-": np.nan, ".": np.nan, "": np.nan}), errors="coerce")
    for c in ["ID", "VARIANT_TYPE", "X_HIN"]:
        df[c] = df[c].fillna(".")
    return df


# ======================================================
# Filtering and measurement
# ======================================================
def apply_filters(df, params):
    """Re-applies the bcftools prefilter, then the gnomAD and VAF filters.

    Note: bcftools evaluates ALT!~"," per allele, so the prefilter keeps
    multiallelic records; they are removed here, not in step 1.
    """
    df = df[df["QUAL"].fillna(-np.inf) >= params.min_qual]
    df = df[~df["ALT"].astype(str).str.contains(",", regex=False)]
    df = df[df["X_HIL"] <= params.max_hmer_length]
    if params.at_only:
        df = df[df["X_HIN"].isin(["A", "T"])]
    gnomad_af = pd.to_numeric(df["GNOMAD_AF"], errors="coerce")
    df = df[(df["ID"] == ".") & ((gnomad_af < params.gnomad_af_threshold) | gnomad_af.isna())]
    df = df[df["VAF"].fillna(np.inf) <= params.max_vaf]
    return df


def _safe_pct(numer, denom):
    return numer / denom * 100 if denom else np.nan


def _safe_ratio(numer, denom):
    return numer / denom if denom else np.nan


def measure_row(df, vcf_name, vcf_path, n_prefiltered, params):
    vt = df["VARIANT_TYPE"].astype(str).to_numpy()
    is_h = vt == "h-indel"
    n_hindel = int(is_h.sum())
    total_indels = n_hindel + int((vt == "non-h-indel").sum())

    ref_len = df["REF"].astype(str).str.len().to_numpy()
    alt_len = df["ALT"].astype(str).str.len().to_numpy()
    one_nt_del = (ref_len - alt_len) == 1
    two_nt_del = (ref_len - alt_len) == 2
    one_nt_ins = (alt_len - ref_len) == 1

    x_hil = df["X_HIL"].to_numpy(dtype=float)
    is_long_h = (x_hil >= params.min_long_hmer_length) & (x_hil <= params.max_hmer_length) & is_h

    n_1nt_del = int((one_nt_del & is_h).sum())
    n_2nt_del = int((two_nt_del & is_h).sum())
    n_1nt_ins = int((one_nt_ins & is_h).sum())
    n_long_hindels = int(is_long_h.sum())
    n_1nt_del_in_long = int((one_nt_del & is_long_h).sum())

    return {
        "vcf": vcf_name,
        "vcf_path": vcf_path,
        "n_records_in_bed": n_prefiltered,
        "n_records_after_filters": len(df),
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
        "1 nt del / 1 nt ins": _safe_ratio(n_1nt_del, n_1nt_ins),
        "2 nt del / 1 nt ins": _safe_ratio(n_2nt_del, n_1nt_ins),
    }


if __name__ == "__main__":
    main()
