#!/usr/bin/env python3
"""Related-element survey across the 159-genome public panel (Supplementary Table S11).

Usage:
    python panel_family_distribution.py --tblastn <marker_vs_panel.tblastn.tsv> \
        --metadata <panel_metadata.tsv> --out <TableS11.tsv>

Inputs
------
* tblastn output (outfmt 6, 12 columns: qseqid sseqid pident length qstart qend
  sstart send evalue bitscore qlen slen) of the four element marker proteins
  (integrase GL001245, large terminase GL001272, portal GL001273, tape-measure
  GL001284) against the panel nucleotide FASTA (BLAST+ 2.17.0, E <= 1e-5).
* panel metadata TSV with columns: accession, organism, strain, assembly_level,
  isolation_source (see results/tables/supplementary/TableS3_public_panel/).
* The panel itself is available from NCBI (accessions in the metadata file);
  genomes are not redistributed in this repository.

A genome is scored as a related-element carrier when its large terminase matches
at >= 60% identity over >= 90% of the query. Result: 19 genomes (18 public plus
the public strain 194 reference itself), large-terminase identity 80.0-98.6%.
"""
import argparse
import csv
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tblastn", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    meta = {r["accession"]: r for r in csv.DictReader(open(args.metadata, encoding="utf-8"), delimiter="\t")}
    rows = [l.rstrip("\n").split("\t") for l in open(args.tblastn)]

    terl, portal = defaultdict(list), defaultdict(list)
    for r in rows:
        q, s = r[0], r[1]
        gcf = s.split("|")[0]
        pident = float(r[2])
        qcov = (int(r[5]) - int(r[4]) + 1) / int(r[10])
        if q == "smbu06GL001272":      # large terminase
            terl[gcf].append((pident, qcov, s))
        elif q == "smbu06GL001273":    # portal
            portal[gcf].append((pident, qcov, s))

    out = []
    for gcf, hits in terl.items():
        best = max(hits, key=lambda x: x[0] * x[1])
        if best[1] >= 0.9 and best[0] >= 60:
            port = max(portal.get(gcf, [(0.0, 0.0, "-")]), key=lambda x: x[0] * x[1])
            m = meta.get(gcf, {})
            out.append(dict(accession=gcf, organism=m.get("organism", "?"), strain=m.get("strain", "?"),
                            assembly_level=m.get("assembly_level", "?"), isolation=m.get("isolation_source", ""),
                            terL_pident=round(best[0], 1), terL_qcov=round(best[1], 2),
                            portal_pident=round(port[0], 1), portal_qcov=round(port[1], 2),
                            contig=best[2]))
    out.sort(key=lambda x: -x["terL_pident"])
    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=list(out[0].keys()))
        w.writeheader()
        for r in out:
            w.writerow(r)
    print(f"{len(out)} family members -> {args.out}")


if __name__ == "__main__":
    main()
