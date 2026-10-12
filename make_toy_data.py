#!/usr/bin/env python3
"""
Generate the synthetic toy data in example/ (install check only; not real data).

Writes:
  example/toy_regions.bed            two small chr1 regions (BED, 0-based)
  example/vcfs/toy_msi.vcf.gz(.tbi)  MSI-like: enriched for 1-nt deletions in long
                                     homopolymers; carries INFO/GNOMAD_AF
  example/vcfs/toy_mss.vcf.gz(.tbi)  MSS-like: no INFO/GNOMAD_AF, so the run
                                     exercises the --gnomad_tsv annotation path
  example/toy_gnomad_af.tsv.gz(.tbi) AFs for toy_mss records

Each VCF also contains records that every filter should remove (non-PASS, low or
missing QUAL, known ID, multiallelic, homopolymer > 15, gnomAD-common, VAF > 0.5,
missing VAF, outside the BED, SNPs).

Requires bgzip and tabix (htslib) on PATH.
"""

import os
import random
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
REGIONS = [("chr1", 100_000, 200_000), ("chr1", 300_000, 400_000)]
OUTSIDE = ("chr1", 250_000, 260_000)

PROFILES = {
    "toy_msi": dict(long_1del=60, short_1del=25, h_1ins=8, h_2del=6, nonh=20, snp=30),
    "toy_mss": dict(long_1del=4, short_1del=10, h_1ins=8, h_2del=3, nonh=25, snp=30),
}
NOISE = ["lowqual_filter", "qual_low", "qual_missing", "known_id", "multiallelic",
         "hil_gt_max", "gnomad_common", "vaf_high", "vaf_missing", "outside_bed"]
N_NOISE_EACH = 3

HEADER = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowQual,Description="Low quality">
##contig=<ID=chr1,length=248956422>
##INFO=<ID=VARIANT_TYPE,Number=1,Type=String,Description="snp, h-indel or non-h-indel">
##INFO=<ID=X_HIL,Number=A,Type=Integer,Description="Homopolymer indel length">
##INFO=<ID=X_HIN,Number=A,Type=String,Description="Homopolymer indel nucleotide">
##INFO=<ID=X_IC,Number=A,Type=String,Description="Indel class">
{gnomad}##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=VAF,Number=A,Type=Float,Description="Variant allele frequency">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}
"""
GNOMAD_HDR = ('##INFO=<ID=GNOMAD_AF,Number=A,Type=Float,'
              'Description="Gnomad Allele Frequency">\n')


def _hbase(rng):
    return rng.choices("ACGT", weights=[0.4, 0.1, 0.1, 0.4])[0]


def make_record(rng, kind):
    """Return dict(ref, alt, info, qual, filt, id, vaf, af, outside)."""
    b = _hbase(rng)
    anchor = "G" if b in "AT" else "A"  # anchor base differs from the repeat base
    rec = dict(qual=rng.randint(20, 60), filt="PASS", id=".",
               vaf=round(rng.uniform(0.05, 0.45), 3),
               af=rng.choice([None, None, round(rng.uniform(1e-5, 5e-3), 5)]),
               outside=False)

    def hindel(ref, alt, hil, ic):
        rec.update(ref=ref, alt=alt,
                   info=f"VARIANT_TYPE=h-indel;X_HIL={hil};X_HIN={b};X_IC={ic}")

    if kind in ("long_1del", "gnomad_common", "vaf_high", "vaf_missing", "known_id",
                "qual_low", "qual_missing", "lowqual_filter", "outside_bed"):
        hindel(anchor + b, anchor, rng.randint(10, 15), "del")
    elif kind == "short_1del":
        hindel(anchor + b, anchor, rng.randint(2, 9), "del")
    elif kind == "h_1ins":
        hindel(anchor, anchor + b, rng.randint(2, 15), "ins")
    elif kind == "h_2del":
        hindel(anchor + b * 2, anchor, rng.randint(4, 15), "del")
    elif kind == "hil_gt_max":
        hindel(anchor + b, anchor, rng.randint(16, 20), "del")
    elif kind == "multiallelic":
        rec.update(ref=anchor + b * 2, alt=f"{anchor + b},{anchor}",
                   info=f"VARIANT_TYPE=h-indel;X_HIL=12,12;X_HIN={b},{b};X_IC=del,del")
        rec["vaf"] = "0.2,0.1"
    elif kind == "nonh":
        seq = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 5)))
        is_del = rng.random() < 0.5
        rec.update(ref=anchor + seq if is_del else anchor,
                   alt=anchor if is_del else anchor + seq,
                   info="VARIANT_TYPE=non-h-indel;X_HIL=0;X_HIN=.;"
                        f"X_IC={'del' if is_del else 'ins'}")
    elif kind == "snp":
        rec.update(ref="A", alt="G", info="VARIANT_TYPE=snp;X_HIL=0;X_HIN=.;X_IC=.")
    else:
        raise ValueError(kind)

    if kind == "lowqual_filter":
        rec["filt"] = "LowQual"
    elif kind == "qual_low":
        rec["qual"] = rng.randint(1, 15)
    elif kind == "qual_missing":
        rec["qual"] = "."
    elif kind == "known_id":
        rec["id"] = f"rs{rng.randint(10**6, 10**8)}"
    elif kind == "gnomad_common":
        rec["af"] = round(rng.uniform(0.02, 0.4), 4)
    elif kind == "vaf_high":
        rec["vaf"] = round(rng.uniform(0.55, 1.0), 3)
    elif kind == "vaf_missing":
        rec["vaf"] = "."
    elif kind == "outside_bed":
        rec["outside"] = True
    return rec


def build_records(profile, seed, noise_each=N_NOISE_EACH):
    rng = random.Random(seed)
    kinds = [k for k, n in profile.items() for _ in range(n)]
    kinds += [k for k in NOISE for _ in range(noise_each)]
    recs = [make_record(rng, k) for k in kinds]

    used = set()
    for r in recs:
        chrom, start, end = OUTSIDE if r["outside"] else rng.choice(REGIONS)
        while True:
            pos = rng.randint(start + 10, end - 10)  # 1-based, well inside the interval
            if pos not in used:
                used.add(pos)
                break
        r["chrom"], r["pos"] = chrom, pos
    return sorted(recs, key=lambda r: r["pos"])


def write_vcf(path, sample, recs, gnomad_in_info):
    lines = [HEADER.format(gnomad=GNOMAD_HDR if gnomad_in_info else "", sample=sample)]
    for r in recs:
        info = r["info"]
        if gnomad_in_info and r["af"] is not None and "," not in r["alt"]:
            info += f";GNOMAD_AF={r['af']}"
        lines.append("\t".join(map(str, [
            r["chrom"], r["pos"], r["id"], r["ref"], r["alt"], r["qual"], r["filt"],
            info, "GT:VAF", f"0/1:{r['vaf']}"])) + "\n")
    plain = path[:-3]
    with open(plain, "w") as fh:
        fh.writelines(lines)
    subprocess.run(["bgzip", "-f", plain], check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", path], check=True)


def write_gnomad_tsv(path, recs):
    rows = sorted((r["chrom"], r["pos"], r["ref"], r["alt"], r["af"])
                  for r in recs if r["af"] is not None and "," not in r["alt"])
    plain = path[:-3]
    with open(plain, "w") as fh:
        for row in rows:
            fh.write("\t".join(map(str, row)) + "\n")
    subprocess.run(["bgzip", "-f", plain], check=True)
    subprocess.run(["tabix", "-f", "-s1", "-b2", "-e2", path], check=True)


def write_bed(path):
    with open(path, "w") as fh:
        for chrom, start, end in REGIONS:
            fh.write(f"{chrom}\t{start}\t{end}\n")


def main():
    os.makedirs(os.path.join(HERE, "vcfs"), exist_ok=True)
    write_bed(os.path.join(HERE, "toy_regions.bed"))

    msi = build_records(PROFILES["toy_msi"], seed=1)
    write_vcf(os.path.join(HERE, "vcfs", "toy_msi.vcf.gz"), "toy_msi", msi,
              gnomad_in_info=True)

    mss = build_records(PROFILES["toy_mss"], seed=2)
    write_vcf(os.path.join(HERE, "vcfs", "toy_mss.vcf.gz"), "toy_mss", mss,
              gnomad_in_info=False)
    write_gnomad_tsv(os.path.join(HERE, "toy_gnomad_af.tsv.gz"), mss)
    print("Toy data written to", HERE)


if __name__ == "__main__":
    main()
