#!/usr/bin/env python3
"""Aggregate tool outputs into Supplementary Table S10 (prophage-tool cross-validation).

This script PARSES the real tool outputs; nothing about the outcome is hardcoded.

Inputs
------
--genomad-root : directory with one subdirectory per geNomad run, e.g.
                 <root>/element/  <root>/smbu06/  <root>/smbu08/  <root>/cc/
                 Each run dir is searched recursively (flattened file names with
                 "__" separators, as published, are fine) for
                 *_virus_summary.tsv, *_marker_classification.tsv and
                 *_nn_classification.tsv. If a run dir is absent -> "input
                 missing"; if present but a needed table is absent -> "run failed
                 or output missing"; if several files match one pattern ->
                 "unparseable (ambiguous match)".
--phispy-root  : directory with one subdirectory per genome, e.g.
                 <root>/smbu06/ (containing phispy_smbu06.log and
                 prophage_coordinates_smbu06.tsv) and <root>/smbu08/.
                 The log is parsed for the ranked candidate table ("Potential
                 prophages") and for the final coordinates file. A genome with no
                 retained call and a logged dropped candidate is reported as
                 "run succeeded - candidate dropped (no retained call)".
--target-chromosome/--target-start/--target-end : the element of interest
                 (1-based inclusive coordinates on that chromosome). Used to
                 compute the actual overlap of each called region with the target;
                 nothing is assumed about which region "should" be called.
--status-file  : optional TSV (tool, input, called, ..., notes) with study-specific
                 tool-status rows (e.g. PHASTER not run; VIBRANT/CheckV not
                 completed). These are facts about a particular study, NOT
                 defaults, so they are supplied explicitly via this file. An
                 example matching the published study is provided as
                 tool_status_example.tsv.

Coordinate convention: all coordinates are 1-based inclusive.
Outcome classes per (tool, run): detected_target_overlap / no_target_overlap /
run_failed / input_missing / unparseable. Negative or failed categories are never
re-labelled as completed results.

Example (published study):
    python make_table_s10.py \
        --genomad-root ../../../results/tables/element_tool_validation/genomad \
        --phispy-root  ../../../results/tables/element_tool_validation/phispy \
        --target-chromosome Chromosome1 --target-start 1255783 --target-end 1294501 \
        --status-file tool_status_example.tsv \
        --out ../../../results/tables/supplementary/TableS10_tool_crossvalidation/TableS10_tool_crossvalidation.tsv
"""
import argparse
import csv
import glob
import os
import re
import sys

COLS = ["tool", "input", "called", "coordinates", "length", "score", "hallmarks",
        "class_", "overlap_with_element", "notes", "raw_output"]


def read_tsv(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def unique_match(root, pattern):
    """Return (path, None) for a unique match; (None, 'missing'/'ambiguous') otherwise."""
    hits = sorted(glob.glob(os.path.join(root, "**", pattern), recursive=True))
    if not hits:
        return None, "missing"
    if len(hits) > 1:
        return None, "ambiguous"
    return hits[0], None


def parse_interval(coord_str):
    """Parse 'start-end' (1-based inclusive) -> (start, end) or None."""
    m = re.match(r"\s*(\d+)\s*-\s*(\d+)\s*$", coord_str or "")
    if not m:
        return None
    return int(m.group(1)), int(m.group(2))


def overlap_bp(a, b):
    """Overlap length of two 1-based inclusive intervals (0 if disjoint)."""
    if a is None or b is None:
        return 0
    return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)


def parse_phispy_log(path):
    """Extract the ranked candidate table from a PhiSpy log.

    Returns dict(candidates=[(contig, start, stop, n_genes, status)], run_failed=bool).
    """
    out = {"candidates": [], "run_failed": False}
    text = open(path, encoding="utf-8", errors="replace").read()
    if re.search(r"Traceback|ModuleNotFoundError|Segmentation fault", text):
        out["run_failed"] = True
    m = re.search(r"Potential prophages \(sorted highest to lowest\)[^\n]*\n", text)
    if m:
        body = text[m.end():]
        for line in body.split("\n"):
            # tolerate both plain stdout and the internal log format
            # ("2026-... INFO     Contig\tStart\t...")
            line = re.sub(r"^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2} \w+\s+", "", line.rstrip("\n"))
            if line.strip() == "" or line.split("\t")[0].strip() == "Contig":
                continue
            parts = line.split("\t")
            if len(parts) >= 5 and parts[1].strip().isdigit():
                out["candidates"].append((parts[0].strip(), int(parts[1]), int(parts[2]),
                                          parts[3].strip(), parts[4].strip()))
            else:
                break
    return out


def parse_phispy_coords(path):
    """Parse prophage_coordinates.tsv -> list of (pp_id, contig, start, stop)."""
    rows = []
    if not path or not os.path.exists(path):
        return rows
    for r in read_tsv(path):
        vals = list(r.values())
        if len(vals) >= 4 and str(vals[2]).strip().isdigit():
            rows.append((vals[0], vals[1], int(vals[2]), int(vals[3])))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genomad-root", required=True)
    ap.add_argument("--phispy-root", required=True)
    ap.add_argument("--target-chromosome", required=True,
                    help="chromosome name of the element of interest (1-based inclusive coordinates)")
    ap.add_argument("--target-start", type=int, required=True)
    ap.add_argument("--target-end", type=int, required=True)
    ap.add_argument("--status-file", default=None,
                    help="optional TSV with study-specific tool-status rows (see tool_status_example.tsv)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    target = (args.target_start, args.target_end)
    rows = []

    # ---------------- geNomad ----------------
    run_desc = {
        "element": "element sequence (attL-attR, 38,865 bp)",
        "smbu06": "complete smbu06 assembly (chromosome + plasmid)",
        "smbu08": "complete smbu08 assembly (chromosome + plasmid + circular component)",
        "cc": "77,437-bp circular component (isolated)",
    }
    for tag, desc in run_desc.items():
        base = os.path.join(args.genomad_root, tag)
        if not os.path.isdir(base):
            rows.append(dict(tool="geNomad v1.12.0", input=desc, called="input_missing",
                             coordinates="", length="", score="", hallmarks="", class_="",
                             overlap_with_element="", notes="run directory absent",
                             raw_output=base))
            continue
        vs, why = unique_match(base, "*_virus_summary.tsv")
        if vs is None:
            rows.append(dict(tool="geNomad v1.12.0", input=desc,
                             called=("unparseable" if why == "ambiguous" else "run_failed"),
                             coordinates="", length="", score="", hallmarks="", class_="",
                             overlap_with_element="", notes=f"virus_summary: {why}",
                             raw_output=base))
        else:
            for r in read_tsv(vs):
                iv = parse_interval(r.get("coordinates", ""))
                seq_level = iv is None  # e.g. whole-sequence virus calls carry no interval
                ov = overlap_bp(iv, target) if iv else 0
                if seq_level:
                    called = "virus (sequence-level; no genomic interval)"
                    ov_txt = "not applicable (whole-sequence input)"
                elif ov > 0:
                    called = "detected_target_overlap"
                    ov_txt = f"{ov} bp ({100.0*ov/(target[1]-target[0]+1):.1f}% of target)"
                else:
                    called = "no_target_overlap"
                    ov_txt = "none"
                rows.append(dict(
                    tool="geNomad v1.12.0 (db v1.9, 2025-02-18; MMseqs2 18.8cc5c)",
                    input=desc, called=called,
                    coordinates=r.get("coordinates", ""), length=r.get("length", ""),
                    score=r.get("virus_score", ""), hallmarks=r.get("n_hallmarks", ""),
                    class_=r.get("taxonomy", ""), overlap_with_element=ov_txt,
                    notes=("provirus call" if r.get("topology") == "Provirus" else "virus (sequence-level)"),
                    raw_output=os.path.relpath(vs).replace("\\", "/")))
        for label, pat in (("geNomad marker-based score", "*_marker_classification.tsv"),
                           ("geNomad NN-based score", "*_nn_classification.tsv")):
            hits = sorted(glob.glob(os.path.join(base, "**", pat), recursive=True))
            # prefer the whole-sequence classification over per-provirus variants
            non_prov = [h for h in hits if "provirus" not in os.path.basename(h)]
            pick = non_prov if len(non_prov) == 1 else (hits if len(hits) == 1 else None)
            if pick is None:
                why = "ambiguous" if hits else "missing"
                rows.append(dict(tool=label, input=desc,
                                 called=("unparseable" if why == "ambiguous" else "run_failed"),
                                 coordinates="", length="", score="", hallmarks="", class_="",
                                 overlap_with_element="", notes=f"classification table: {why}",
                                 raw_output=base))
                continue
            f = pick[0]
            for r in read_tsv(f):
                rows.append(dict(tool=label, input=desc, called="score detail",
                                 coordinates=r.get("seq_name", ""), length="",
                                 score=(f"virus={r.get('virus_score')}; chromosome={r.get('chromosome_score')}; "
                                        f"plasmid={r.get('plasmid_score')}"),
                                 hallmarks="", class_="", overlap_with_element="", notes="",
                                 raw_output=os.path.relpath(f).replace("\\", "/")))

    # ---------------- PhiSpy ----------------
    genome_desc = {
        "smbu06": "complete smbu06 assembly, GenBank (provider annotation)",
        "smbu08": "complete smbu08 assembly, GenBank (provider annotation)",
    }
    for tag, desc in genome_desc.items():
        # layout A: <root>/<tag>/...   layout B (as published): flat files in <root>
        gdir = os.path.join(args.phispy_root, tag)
        if os.path.isdir(gdir):
            log_hits = sorted(glob.glob(os.path.join(gdir, "*.log")))
            coords_pat = os.path.join(gdir, "*coordinates*.tsv")
        else:
            gdir = args.phispy_root
            log_hits = sorted(glob.glob(os.path.join(gdir, f"*{tag}*.log")))
            coords_pat = os.path.join(gdir, f"*coordinates*{tag}*.tsv")
        if not log_hits:
            rows.append(dict(tool="PhiSpy v5.0.10 (Akhter et al. 2012)", input=desc,
                             called="input_missing", coordinates="", length="", score="",
                             hallmarks="", class_="", overlap_with_element="",
                             notes="log directory absent", raw_output=gdir))
            continue
        log = log_hits[0]
        parsed = parse_phispy_log(log)
        if parsed["run_failed"]:
            rows.append(dict(tool="PhiSpy v5.0.10 (Akhter et al. 2012)", input=desc,
                             called="run_failed", coordinates="", length="", score="",
                             hallmarks="", class_="", overlap_with_element="",
                             notes="error markers in log", raw_output=os.path.relpath(log).replace("\\", "/")))
            continue
        # target-overlapping candidates (any status)
        tgt_cands = [c for c in parsed["candidates"]
                     if c[0] == args.target_chromosome and overlap_bp((c[1], c[2]), target) > 0]
        coords_files = sorted(glob.glob(coords_pat))
        kept = parse_phispy_coords(coords_files[0]) if coords_files else []
        kept_tgt = [k for k in kept if k[1] == args.target_chromosome and overlap_bp((k[2], k[3]), target) > 0]
        if kept_tgt:
            for k in kept_tgt:
                ov = overlap_bp((k[2], k[3]), target)
                rows.append(dict(tool="PhiSpy v5.0.10", input=desc, called="detected_target_overlap",
                                 coordinates=f"{k[2]}-{k[3]}", length=str(k[3] - k[2] + 1),
                                 score="", hallmarks="", class_="",
                                 overlap_with_element=f"{ov} bp", notes=f"retained call {k[0]}",
                                 raw_output=os.path.relpath(log).replace("\\", "/")))
        elif tgt_cands:
            det = "; ".join(f"{c[0]}:{c[1]}-{c[2]} ({c[3]} genes; {c[4]})" for c in tgt_cands)
            rows.append(dict(tool="PhiSpy v5.0.10", input=desc,
                             called="no_target_overlap (candidate dropped; no retained call)",
                             coordinates=det, length="", score="",
                             hallmarks="", class_="",
                             overlap_with_element="candidate overlaps the target region",
                             notes=("the ranked candidate overlapping the target was not retained; "
                                    "see the drop message in the coordinates field and the run log"),
                             raw_output=os.path.relpath(log).replace("\\", "/")))
        else:
            rows.append(dict(tool="PhiSpy v5.0.10", input=desc,
                             called=("no_target_overlap (run succeeded; no candidate at/near the target)"
                                     if parsed["candidates"] else "no_target_overlap (run succeeded; no candidates logged)"),
                             coordinates="", length="", score="", hallmarks="", class_="",
                             overlap_with_element="",
                             notes="no candidate region overlapping the target interval",
                             raw_output=os.path.relpath(log).replace("\\", "/")))
        # retained calls elsewhere (reported, not converted into anything)
        other_kept = [k for k in kept if k not in kept_tgt]
        for k in other_kept:
            rows.append(dict(tool="PhiSpy v5.0.10", input=desc,
                             called="no_target_overlap (retained call elsewhere)",
                             coordinates=f"{k[1]} {k[2]}-{k[3]}", length=str(k[3] - k[2] + 1),
                             score="", hallmarks="", class_="", overlap_with_element="none",
                             notes=f"retained call {k[0]} is outside the target interval",
                             raw_output=os.path.relpath(log).replace("\\", "/")))

    # ---------------- study-specific status rows (explicit config, not defaults) ----------------
    if args.status_file:
        if not os.path.exists(args.status_file):
            sys.exit(f"status file not found: {args.status_file}")
        for r in read_tsv(args.status_file):
            rows.append({c: r.get(c, "") for c in COLS})

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=COLS)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"{len(rows)} rows -> {args.out}")


if __name__ == "__main__":
    main()
