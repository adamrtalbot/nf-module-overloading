#!/usr/bin/env python3
"""Compare primary alignments in several SAM files against a baseline.

usage: compare_sams.py BASELINE_ID id1.sam id2.sam ...  (id is the file name without .sam)
"""
import sys
from pathlib import Path

FIELDS = {"flag": 1, "chr": 2, "pos": 3, "mapq": 4, "cigar": 5}


def primary_alignments(path):
    """Primary records keyed by (read name, mate)."""
    records = {}
    for line in open(path):
        if line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        flag = int(f[1])
        if flag & 0x900:  # secondary or supplementary
            continue
        mate = 2 if flag & 0x80 else 1
        records[(f[0], mate)] = f
    return records


baseline, files = sys.argv[1], sys.argv[2:]
parsed = {Path(p).stem: primary_alignments(p) for p in files}
base = parsed[baseline]

cols = ["records", "mapped", "chr=", "pos=", "cigar=", "mapq="]
print(f"Alignments compared to {baseline} (primary records matched by read and mate)")
print(f"{'aligner':<10}" + "".join(f"{c:>9}" for c in cols))
for name in sorted(parsed):
    aln = parsed[name]
    mapped = sum(1 for f in aln.values() if not int(f[1]) & 0x4)
    same = {
        k: sum(1 for key, b in base.items() if key in aln and aln[key][i] == b[i])
        for k, i in FIELDS.items()
    }
    vals = [len(aln), mapped, same["chr"], same["pos"], same["cigar"], same["mapq"]]
    print(f"{name:<10}" + "".join(f"{v:>9}" for v in vals))
