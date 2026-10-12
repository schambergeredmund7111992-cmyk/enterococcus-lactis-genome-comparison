#!/usr/bin/env python3
"""Aggregate geNomad and PhiSpy outputs into Supplementary Table S10
(prophage-tool cross-validation).

Usage:
    python make_table_s10.py --genomad-root <dir with per-run geNomad outputs> \
        --phispy-root <dir with phispy logs / coordinates> \
        --out results/tables/supplementary/TableS10_tool_crossvalidation/TableS10_tool_crossvalidation.tsv

The expected layout under --genomad-root mirrors the published runs:
    <root>/element/  <root>/smbu06/  <root>/smbu08/  <root>/cc/
each containing *_virus_summary.tsv, *_provirus.tsv, *_marker_classification.tsv,
*_nn_classification.tsv (flattened file names are accepted).

Tool status rows (PHASTER not run; VIBRANT/CheckV not completed; provider tool
not disclosed) are written as-is and must not be reported as negative results.
"""
import argparse
import csv
import glob
import os


def read_tsv(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def find_one(root, pattern):
    hits = glob.glob(os.path.join(root, "**", pattern), recursive=True)
    return hits[0] if hits else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genomad-root", required=True)
    ap.add_argument("--phispy-root", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rows = []
    run_desc = {
        "element": "element sequence (attL-attR, 38,865 bp)",
        "smbu06": "complete smbu06 assembly (chromosome + plasmid)",
        "smbu08": "complete smbu08 assembly (chromosome + plasmid + circular component)",
        "cc": "77,437-bp circular component (isolated)",
    }
    for tag, desc in run_desc.items():
        base = os.path.join(args.genomad_root, tag)
        if not os.path.isdir(base):
            rows.append(dict(tool="geNomad v1.12.0", input=desc, called="output missing",
                             coordinates="", length="", score="", hallmarks="", class_="",
                             overlap_with_element="", notes="", raw_output=base))
            continue
        vs = find_one(base, "*_virus_summary.tsv")
        if vs:
            for r in read_tsv(vs):
                rows.append(dict(
                    tool="geNomad v1.12.0 (db v1.9, 2025-02-18; MMseqs2 18.8cc5c)",
                    input=desc,
                    called=("provirus call" if r.get("topology") == "Provirus" else "virus (sequence-level)"),
                    coordinates=r.get("coordinates", ""), length=r.get("length", ""),
                    score=r.get("virus_score", ""), hallmarks=r.get("n_hallmarks", ""),
                    class_=r.get("taxonomy", ""), overlap_with_element="", notes="",
                    raw_output=os.path.relpath(vs).replace("\\", "/")))
        for label, pat in (("geNomad marker-based score", "*_marker_classification.tsv"),
                           ("geNomad NN-based score", "*_nn_classification.tsv")):
            f = find_one(base, pat)
            if f:
                for r in read_tsv(f):
                    rows.append(dict(tool=label, input=desc, called="score detail",
                                     coordinates=r.get("seq_name", ""), length="",
                                     score=(f"virus={r.get('virus_score')}; chromosome={r.get('chromosome_score')}; "
                                            f"plasmid={r.get('plasmid_score')}"),
                                     hallmarks="", class_="", overlap_with_element="", notes="",
                                     raw_output=os.path.relpath(f).replace("\\", "/")))

    # PhiSpy summary rows (edit to match the actual run)
    ph = args.phispy_root
    rows.append(dict(tool="PhiSpy v5.0.10 (Akhter et al. 2012)",
                     input="complete smbu06 assembly, GenBank (provider annotation)",
                     called="no retained prophage call (top-ranked candidate dropped)",
                     coordinates="top-ranked candidate 1,264,149-1,279,081 (24 genes, overlaps the element core); dropped with the message 'No genes were identified as phage genes'",
                     length="", score="(no numeric score; candidates are ranked, then filtered by the gene-level phage-gene classification)",
                     hallmarks="", class_="",
                     overlap_with_element="candidate overlaps the element core (1,255,783-1,294,501)",
                     notes="PhiSpy evaluation keeps a candidate only if at least one of its genes is classified as a phage gene (default threshold 1; evaluation.py v5.0.10)",
                     raw_output=os.path.relpath(ph).replace("\\", "/")))
    rows.append(dict(tool="PhiSpy v5.0.10",
                     input="complete smbu08 assembly, GenBank (provider annotation)",
                     called="1 retained prophage call (NOT the study element)",
                     coordinates="Plasmid1 110,154-128,787 (IS/recombinase-rich region of the shared plasmid)",
                     length="18,634", score="", hallmarks="", class_="", overlap_with_element="no",
                     notes="circular-component regions were high-ranked candidates, dropped for the same reason as in smbu06; no candidate at the element position on the chromosome",
                     raw_output=os.path.relpath(ph).replace("\\", "/")))
    rows.append(dict(tool="provider prophage tool (name/version not disclosed in delivery)",
                     input="smbu06/smbu08 assemblies", called="no prophage called (empty table)",
                     coordinates="", length="", score="", hallmarks="", class_="",
                     overlap_with_element="none",
                     notes="tool identity, version, database and thresholds are not documented in the delivery",
                     raw_output="provider delivery (not included in this repository)"))
    rows.append(dict(tool="PHASTER (web)", input="submission package prepared",
                     called="NOT RUN in this revision (web-only service; authors submit)",
                     coordinates="", length="", score="", hallmarks="", class_="",
                     overlap_with_element="", notes="submission package accompanies the manuscript revision (contains assembly sequences; not in this repository)",
                     raw_output=""))
    rows.append(dict(tool="VIBRANT v1.2.1", input="installed; database setup blocked",
                     called="NOT COMPLETED (Pfam/KEGG sources served over retired FTP endpoints; remaining download measured 0-0.17 MB/s from this network)",
                     coordinates="", length="", score="", hallmarks="", class_="",
                     overlap_with_element="", notes="", raw_output="see docs/element_tool_validation_methods.md"))
    rows.append(dict(tool="CheckV v1.1.1", input="installed",
                     called="NOT COMPLETED (database download measured 0.2-5.5 KB/s from the official sources)",
                     coordinates="", length="", score="", hallmarks="", class_="",
                     overlap_with_element="", notes="CheckV is a quality/completeness assessment, not an identification tool",
                     raw_output="see docs/element_tool_validation_methods.md"))

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"{len(rows)} rows -> {args.out}")


if __name__ == "__main__":
    main()
