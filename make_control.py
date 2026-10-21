#!/usr/bin/env python3
import sys
import glob
import os

SEX_FILE   = "./"
READS_DIR  = "/01_reads"
OUT_FILE   = "./control/tsv"
DELIM = " "   # single space; use "  " for double space, "\t" for tab

# --- Index fastq files in 01_reads by (numeric WES id, barcode) ---
# e.g. 11720_EX_WES_289_216508227_S1_R1_001.fastq.gz
#      -> key (289, "216508227") -> (prefix "11720", R1 filename, R2 filename)
samples = {}
r1_files = glob.glob(os.path.join(READS_DIR, "*_R1_001.fastq.gz"))

for r1_path in r1_files:
    r1_name = os.path.basename(r1_path)
    r2_name = r1_name.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz")
    r2_path = os.path.join(READS_DIR, r2_name)

    if not os.path.exists(r2_path):
        print(f"WARNING: missing R2 for {r1_name}", file=sys.stderr)
        continue

    stem = r1_name[:-len("_R1_001.fastq.gz")]   # 11720_EX_WES_289_216508227_S1
    p = stem.split("_")                          # ['11720','EX','WES','289','216508227','S1']
    prefix, num, barcode = p[0], p[3], p[4]

    key = (int(num), barcode)
    if key in samples:
        print(f"WARNING: duplicate fastq entry for key {key}", file=sys.stderr)
    samples[key] = (prefix, r1_name, r2_name)

# --- Parse sex info and build output ---
females, males = [], []
family_lines, fastq_lines = [], []

with open(SEX_FILE) as f:
    next(f)  # skip header line "ID Sex"
    for line in f:
        if not line.strip():
            continue
        sample_id, sex = line.split()[:2]
        p = sample_id.split("_")          # EX_WES_00865_216510662
        padded_num, barcode = p[2], p[3]

        key = (int(padded_num), barcode)
        if key not in samples:
            print(f"WARNING: no fastq match in {READS_DIR} for {sample_id}", file=sys.stderr)
            continue

        prefix, r1, r2 = samples[key]
        sample = f"{prefix}_EX_WES_{padded_num}"   # e.g. 11720_EX_WES_00865

        if sex.upper() == "FEMALE":
            females.append(sample)
        elif sex.upper() == "MALE":
            males.append(sample)
        else:
            print(f"WARNING: unexpected sex '{sex}' for {sample_id}", file=sys.stderr)
            continue

        family_lines.append(DELIM.join(["FAMILY", sample]))
        fastq_lines.append(DELIM.join(["FASTQ", sample, "ILLUMINA", barcode, r1, r2]))

with open(OUT_FILE, "w") as out:
    if females:
        out.write(DELIM.join(["FEMALE"] + females) + "\n")
    if males:
        out.write(DELIM.join(["MALE"] + males) + "\n")
    for line in family_lines:
        out.write(line + "\n")
    for line in fastq_lines:
        out.write(line + "\n")

print(f"Wrote {len(family_lines)} samples "
      f"({len(females)} female, {len(males)} male) to {OUT_FILE}")
