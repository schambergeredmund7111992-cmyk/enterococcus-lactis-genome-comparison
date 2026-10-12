#!/usr/bin/env python3
"""CDS-set assignment audit for the element region (Supplementary Table S13).

Usage:
    python make_table_s13_cds_audit.py --element-genes <prophage_genes.tsv> --out <TableS13.tsv>

Input: the provider-annotated element gene table (one row per CDS with columns
locus_tag, start, end, strand, length_aa, product; coordinates 1-based on the
smbu06 chromosome). The delivered table used here is derived from the provider
GenBank record (not redistributed; see input_manifest.tsv).

Regions (1-based inclusive, verified against the delivered FASTA):
    internal sequence     1,255,929-1,294,501  (38,573 bp)
    complete attL-attR    1,255,783-1,294,647  (38,865 bp)

Result (published): 50 of the 51 provider-annotated CDS lie fully inside the
internal sequence; the 51st (smbu06GL001295, 114 bp) lies in the attR region and
crosses the 3' end of the complete attL-attR element. 'the 51 element proteins'
denotes this provider-associated comparison set, not the number of complete CDS
inside the internal sequence.
"""
import argparse
import csv

INTERNAL = (1255929, 1294501)
FULL = (1255783, 1294647)


def rel(s, e, region):
    rs, re_ = region
    if s >= rs and e <= re_:
        return "fully inside"
    if e < rs or s > re_:
        return "no overlap"
    return "overlaps boundary"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--element-genes", required=True, help="prophage_genes.tsv (provider annotation)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    out = []
    for r in csv.DictReader(open(args.element_genes, encoding="utf-8"), delimiter="\t"):
        s, e = int(r["start"]), int(r["end"])
        ri, rf = rel(s, e, INTERNAL), rel(s, e, FULL)
        if ri == "fully inside":
            cls = "A. fully inside internal sequence (and within the attL-attR span)"
        elif rf == "overlaps boundary":
            cls = "C. crosses the attR/attL-attR end boundary (outside the internal sequence)"
        else:
            cls = "D. outside both regions"
        out.append(dict(gene=r["locus_tag"], start=s, end=e, strand=r["strand"],
                        length_bp=e - s + 1, length_aa=r["length_aa"],
                        vs_internal=ri, vs_complete_attL_attR=rf, category=cls,
                        included_in_51_protein_set=("yes (boundary member of the provider-associated set)" if cls.startswith("C") else "yes"),
                        product=r["product"][:60]))
    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=list(out[0].keys()))
        w.writeheader()
        for r in out:
            w.writerow(r)
    from collections import Counter
    print(Counter(r["category"][0] for r in out), "->", args.out)


if __name__ == "__main__":
    main()
