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
                 *_nn_classification.tsv. Outcomes are classed as
                 detected_target_overlap / no_target_overlap /
                 virus (sequence-level; no genomic interval) / run_failed /
                 input_missing / unparseable. An existing-but-empty summary is
                 reported as its own state, never as a silent success.
--phispy-root  : directory with one subdirectory per genome, e.g.
                 <root>/smbu06/ (log + prophage_coordinates tsv) and
                 <root>/smbu08/. Flat layouts (files named *<tag>*) are also
                 accepted. Several matching logs/coordinate files are reported
                 as ambiguous rather than silently picking the first.
                 Log parsing requires PhiSpy's completion marker
                 ("Done!!! Output is in ...") before any outcome is called a
                 success; truncated logs are reported as incomplete.
                 Coordinate files are parsed in PhiSpy's ACTUAL format
                 (v5 writers.py): no header, first four columns
                 = pp-id / contig / start / stop, 1-based inclusive. A header
                 line is detected explicitly (if present) and skipped only in
                 that case; malformed lines raise a descriptive error naming
                 the file and line number and are never silently dropped.
--target-chromosome/--target-start/--target-end : the element of interest
                 (1-based inclusive coordinates on that chromosome). Used to
                 compute the actual overlap of each called region with the
                 target; nothing is assumed about which region "should" be
                 called.
--status-file  : optional TSV (tool, input, called, ..., notes) with study-specific
                 tool-status rows (e.g. PHASTER not run; VIBRANT/CheckV not
                 completed). These are facts about a particular study, NOT
                 defaults, so they are supplied explicitly via this file. An
                 example matching the published study is provided as
                 tool_status_example.tsv.

Coordinate convention: all coordinates are 1-based inclusive. Outcome classes are
never mixed: success-with-negative-result, dropped candidate, failed run,
incomplete run, missing input and unparseable output are distinct states.

Example (published study):
    python make_table_s10.py \
        --genomad-root ../../results/tables/element_tool_validation/genomad \
        --phispy-root  ../../results/tables/element_tool_validation/phispy \
        --target-chromosome Chromosome1 --target-start 1255783 --target-end 1294501 \
        --status-file tool_status_example.tsv \
        --out ../../results/tables/supplementary/TableS10_tool_crossvalidation/TableS10_tool_crossvalidation.tsv
"""
import argparse
import csv
import glob
import os
import re
import sys

COLS = ["tool", "input", "called", "coordinates", "length", "score", "hallmarks",
        "class_", "overlap_with_element", "notes", "raw_output"]

PHISPY_DONE_RE = re.compile(r"Done!+\s*Output is in", re.I)
PHISPY_FAIL_RE = re.compile(r"Traceback|ModuleNotFoundError|Segmentation fault|"
                            r"unrecognized arguments|No such file or directory")


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
    """Parse a PhiSpy log.

    Returns dict(candidates=[(contig, start, stop, n_genes, status)],
                run_failed=bool, completed=bool).
    'completed' requires PhiSpy's own completion marker ("Done!!! Output is
    in ..."); the mere absence of error keywords is NOT treated as success.
    """
    out = {"candidates": [], "run_failed": False, "completed": False}
    text = open(path, encoding="utf-8", errors="replace").read()
    if PHISPY_FAIL_RE.search(text):
        out["run_failed"] = True
    if PHISPY_DONE_RE.search(text):
        out["completed"] = True
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
    """Parse a PhiSpy prophage_coordinates.tsv file.

    Returns dict(rows=[(pp_id, contig, start, stop)], state=..., detail=...)
    with state in {'ok', 'empty', 'missing', 'malformed'}.

    PhiSpy v5 writes this file WITHOUT a header (writers.py: first four columns
    are pp-id, contig, start, stop). A header line is detected explicitly when
    present and skipped only in that case. Malformed lines are reported with
    file name and line number and are never silently dropped.
    """
    if not path:
        return {"rows": [], "state": "missing", "detail": "no coordinates file located"}
    if not os.path.exists(path):
        return {"rows": [], "state": "missing", "detail": f"file not found: {path}"}
    lines = [l.rstrip("\r\n") for l in open(path, encoding="utf-8", errors="replace")]
    lines = [l for l in lines if l.strip() != ""]
    if not lines:
        return {"rows": [], "state": "empty", "detail": "file is empty (no retained calls)"}

    def is_call(cols):
        return (len(cols) >= 4
                and cols[2].strip().isdigit() and cols[3].strip().isdigit())

    start_idx = 0
    first = lines[0].split("\t")
    if not is_call(first):
        header_words = {"pp", "contig", "start", "stop", "atts"}
        if any(t.strip().lower() in header_words for t in first[:5]):
            start_idx = 1  # explicit header detected
        else:
            return {"rows": [], "state": "malformed",
                    "detail": f"{path}: line 1 is neither a header nor a valid call: {first[:4]!r}"}
    rows = []
    for n, line in enumerate(lines[start_idx:], start=start_idx + 1):
        parts = line.split("\t")
        if not is_call(parts):
            return {"rows": [], "state": "malformed",
                    "detail": f"{path}: malformed call at line {n}: {line[:80]!r}"}
        s, e = int(parts[2]), int(parts[3])
        if s < 1 or e < s:
            return {"rows": [], "state": "malformed",
                    "detail": f"{path}: invalid 1-based coordinates at line {n}: {s}-{e}"}
        rows.append((parts[0].strip(), parts[1].strip(), s, e))
    return {"rows": rows, "state": "ok", "detail": f"{len(rows)} retained call(s)"}


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
            vs_rows = read_tsv(vs)
            if not vs_rows:
                rows.append(dict(tool="geNomad v1.12.0", input=desc,
                                 called="run_failed or empty output",
                                 coordinates="", length="", score="", hallmarks="", class_="",
                                 overlap_with_element="",
                                 notes="virus_summary exists but contains no data rows",
                                 raw_output=os.path.relpath(vs).replace("\\", "/")))
            for r in vs_rows:
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
            cls_rows = read_tsv(f)
            if not cls_rows:
                rows.append(dict(tool=label, input=desc, called="run_failed or empty output",
                                 coordinates="", length="", score="", hallmarks="", class_="",
                                 overlap_with_element="",
                                 notes="classification table exists but contains no data rows",
                                 raw_output=os.path.relpath(f).replace("\\", "/")))
            for r in cls_rows:
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
            coords_hits = sorted(glob.glob(os.path.join(gdir, "*coordinates*.tsv")))
        else:
            gdir = args.phispy_root
            log_hits = sorted(glob.glob(os.path.join(gdir, f"*{tag}*.log")))
            coords_hits = sorted(glob.glob(os.path.join(gdir, f"*coordinates*{tag}*.tsv")))

        if not log_hits:
            rows.append(dict(tool="PhiSpy v5.0.10 (Akhter et al. 2012)", input=desc,
                             called="input_missing", coordinates="", length="", score="",
                             hallmarks="", class_="", overlap_with_element="",
                             notes="no log file found", raw_output=gdir))
            continue
        if len(log_hits) > 1:
            rows.append(dict(tool="PhiSpy v5.0.10 (Akhter et al. 2012)", input=desc,
                             called="unparseable", coordinates="", length="", score="",
                             hallmarks="", class_="", overlap_with_element="",
                             notes=f"multiple log files match ({len(log_hits)}); ambiguous",
                             raw_output="; ".join(os.path.relpath(h).replace("\\", "/") for h in log_hits[:4])))
            continue
        if len(coords_hits) > 1:
            rows.append(dict(tool="PhiSpy v5.0.10 (Akhter et al. 2012)", input=desc,
                             called="unparseable", coordinates="", length="", score="",
                             hallmarks="", class_="", overlap_with_element="",
                             notes=f"multiple coordinates files match ({len(coords_hits)}); ambiguous",
                             raw_output="; ".join(os.path.relpath(h).replace("\\", "/") for h in coords_hits[:4])))
            continue
        log = log_hits[0]
        rawo = os.path.relpath(log).replace("\\", "/")
        parsed = parse_phispy_log(log)

        if parsed["run_failed"]:
            rows.append(dict(tool="PhiSpy v5.0.10 (Akhter et al. 2012)", input=desc,
                             called="run_failed", coordinates="", length="", score="",
                             hallmarks="", class_="", overlap_with_element="",
                             notes="error markers in log", raw_output=rawo))
            continue

        coord_path = coords_hits[0] if coords_hits else None
        coords = parse_phispy_coords(coord_path)

        # --- run-state gates: completion marker and output consistency ---
        if not parsed["completed"]:
            note = ("log lacks PhiSpy's completion marker ('Done!!! Output is in ...'); "
                    "outcome cannot be confirmed (truncated or incomplete run)")
            if coords["state"] == "ok" and coords["rows"]:
                note = ("coordinates file contains calls but the log lacks PhiSpy's completion "
                        "marker; log and output are inconsistent - outcome cannot be confirmed")
            rows.append(dict(tool="PhiSpy v5.0.10 (Akhter et al. 2012)", input=desc,
                             called="incomplete (completion not confirmed)",
                             coordinates="", length="", score="", hallmarks="", class_="",
                             overlap_with_element="", notes=note,
                             raw_output=rawo + (f"; {os.path.relpath(coord_path).replace(chr(92), '/')}" if coord_path else "")))
            continue
        if coords["state"] == "missing":
            cand_note = ("candidate region(s) are listed in the log but the final coordinates "
                         "file is missing; retention status cannot be inferred"
                         if parsed["candidates"] else "no coordinates file and no candidates in log")
            rows.append(dict(tool="PhiSpy v5.0.10 (Akhter et al. 2012)", input=desc,
                             called="incomplete (coordinates file missing)",
                             coordinates="", length="", score="", hallmarks="", class_="",
                             overlap_with_element="", notes=cand_note, raw_output=rawo))
            continue
        if coords["state"] == "malformed":
            rows.append(dict(tool="PhiSpy v5.0.10 (Akhter et al. 2012)", input=desc,
                             called="unparseable", coordinates="", length="", score="",
                             hallmarks="", class_="", overlap_with_element="",
                             notes=coords["detail"], raw_output=rawo))
            continue

        # --- completed run with usable (possibly empty) coordinates file ---
        coords_rel = os.path.relpath(coord_path).replace("\\", "/")
        kept = coords["rows"]
        tgt_cands = [c for c in parsed["candidates"]
                     if c[0] == args.target_chromosome and overlap_bp((c[1], c[2]), target) > 0]
        kept_tgt = [k for k in kept if k[1] == args.target_chromosome
                    and overlap_bp((k[2], k[3]), target) > 0]
        if kept_tgt:
            for k in kept_tgt:
                ov = overlap_bp((k[2], k[3]), target)
                rows.append(dict(tool="PhiSpy v5.0.10", input=desc,
                                 called="detected_target_overlap",
                                 coordinates=f"{k[2]}-{k[3]}", length=str(k[3] - k[2] + 1),
                                 score="", hallmarks="", class_="",
                                 overlap_with_element=f"{ov} bp", notes=f"retained call {k[0]}",
                                 raw_output=rawo + "; " + coords_rel))
        elif tgt_cands:
            det = "; ".join(f"{c[0]}:{c[1]}-{c[2]} ({c[3]} genes; {c[4]})" for c in tgt_cands)
            rows.append(dict(tool="PhiSpy v5.0.10", input=desc,
                             called="no_target_overlap (candidate dropped; no retained call)",
                             coordinates=det, length="", score="", hallmarks="", class_="",
                             overlap_with_element="candidate overlaps the target region",
                             notes=("the ranked candidate overlapping the target was not retained; "
                                    "the coordinates file contains no retained call for it"),
                             raw_output=rawo + "; " + coords_rel))
        else:
            rows.append(dict(tool="PhiSpy v5.0.10", input=desc,
                             called=("no_target_overlap (run completed; no candidate at/near the target)"
                                     if parsed["candidates"] else
                                     "no_target_overlap (run completed; no candidates logged)"),
                             coordinates="", length="", score="", hallmarks="", class_="",
                             overlap_with_element="",
                             notes="no candidate region overlapping the target interval",
                             raw_output=rawo + "; " + coords_rel))
        other_kept = [k for k in kept if k not in kept_tgt]
        for k in other_kept:
            rows.append(dict(tool="PhiSpy v5.0.10", input=desc,
                             called="no_target_overlap (retained call elsewhere)",
                             coordinates=f"{k[1]} {k[2]}-{k[3]}", length=str(k[3] - k[2] + 1),
                             score="", hallmarks="", class_="", overlap_with_element="none",
                             notes=f"retained call {k[0]} is outside the target interval",
                             raw_output=rawo + "; " + coords_rel))

    # ---------------- study-specific status rows (explicit config, not defaults) ----------------
    if args.status_file:
        if not os.path.exists(args.status_file):
            sys.exit(f"status file not found: {args.status_file}")
        for r in read_tsv(args.status_file):
            rows.append({c: r.get(c, "") for c in COLS})

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=COLS, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"{len(rows)} rows -> {args.out}")


if __name__ == "__main__":
    main()
