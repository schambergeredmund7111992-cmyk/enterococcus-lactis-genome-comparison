#!/usr/bin/env python3
"""Self-test for make_table_s10.py using small synthetic fixtures.

Covers (each in an isolated temporary root; no real study data are touched):

geNomad aggregation
  A. target detected      - provirus interval overlapping the target
  D. no target overlap    - call disjoint from the target; sequence-level input flagged
  C. input missing        - expected run directory absent
  E. empty summary        - summary exists but has no data rows (not a silent success)

PhiSpy coordinate parsing (actual headerless format; header tolerated if present)
  F. headerless single call over the target      -> detected_target_overlap
  G. headerless single call outside the target   -> retained call elsewhere
  H. multiple calls (target + elsewhere)         -> both reported
  I. file WITH a header                          -> header detected, rows kept
  J. empty coordinates file + completed run      -> succeeded negative (no retained calls)
  K. malformed coordinates line                  -> unparseable with file/line in notes

PhiSpy run-state gates
  L. truncated log (no completion marker)        -> incomplete, NOT a success-negative
  M. completed run, no candidates                -> succeeded negative
  N. completed run, candidates present, coordinates file MISSING
                                                 -> incomplete; retention not inferred
  O. explicit failure markers                    -> run_failed
  P. candidates + retained rows but NO completion marker (log/output inconsistent)
                                                 -> incomplete; not confirmed

Usage:  python selftest_make_table_s10.py
"""
import csv
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "make_table_s10.py")
TARGET = ("Chromosome1", 1000, 3000)
DONE = "Done!!! Output is in phispy_out and the log is in phispy_out\\phispy.log\n"
CAND_HDR = "Potential prophages (sorted highest to lowest)\nContig\tStart\tStop\tNumber of potential genes\tStatus\n"
CAND_TGT = "Chromosome1\t1200\t2800\t24\tDropped. No genes were identified as phage genes\n"

GSM_HDR = ("seq_name\tlength\ttopology\tcoordinates\tn_genes\tgenetic_code\tvirus_score\t"
           "fdr\tn_hallmarks\tmarker_enrichment\ttaxonomy\n")


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="\n") as f:
        f.write(text)


def run(root):
    out = os.path.join(root, "s10.tsv")
    cmd = [sys.executable, SCRIPT,
           "--genomad-root", os.path.join(root, "genomad"),
           "--phispy-root", os.path.join(root, "phispy"),
           "--target-chromosome", TARGET[0],
           "--target-start", str(TARGET[1]), "--target-end", str(TARGET[2]),
           "--out", out]
    r = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace")
    if r.returncode != 0:
        print("STDERR:", r.stderr)
        raise SystemExit("script failed")
    with open(out, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def fresh():
    return tempfile.mkdtemp(prefix="s10t_")


def phispy_tag(rows, cands=True, done=True, coord_name="prophage_coordinates_smbu06.tsv",
               no_cand_table=False):
    """Build (log_text, coords_text) for one tag."""
    if no_cand_table:
        log = ""
    else:
        log = CAND_HDR + (CAND_TGT if cands else "Plasmid1\t200\t400\t2\tDropped. Region too small (Not enough genes)\n")
    if done:
        log += DONE
    return log, rows, coord_name


def phispy_case(expect_calls, expect_status_sub, rows, cands=True, done=True,
                write_coords=True, coord_name=None, extra_assert=None,
                no_cand_table=False):
    root = fresh()
    try:
        log, rows_txt, cname = phispy_tag(rows, cands=cands, done=done,
                                          no_cand_table=no_cand_table)
        cname = coord_name or cname
        write(os.path.join(root, "phispy", "phispy_smbu06.log"), log)
        if write_coords:
            write(os.path.join(root, "phispy", cname), rows_txt)
        table = run(root)
        hits = [r for r in table if r["tool"].startswith("PhiSpy") and "smbu06" in r["input"]]
        assert hits, f"no PhiSpy smbu06 rows: {table}"
        calls = [h["called"] for h in hits]
        for sub in expect_status_sub:
            assert any(sub in c for c in calls), f"expected {sub!r} in {calls}"
        if extra_assert:
            extra_assert(hits)
        return hits
    finally:
        shutil.rmtree(root, ignore_errors=True)


def main():
    print("self-test: make_table_s10.py (v3)")

    # ---------------- geNomad ----------------
    root = fresh()
    try:
        table = run(root)
        miss = [r for r in table if r["called"] == "input_missing"]
        assert miss and all(r["coordinates"] == "" and r["score"] == "" for r in miss)
        print("C (geNomad input missing)          OK")
    finally:
        shutil.rmtree(root, ignore_errors=True)

    root = fresh()
    try:
        write(os.path.join(root, "genomad", "smbu06", "run_virus_summary.tsv"),
              GSM_HDR + "Chr|provirus_950_3100\t2150\tProvirus\t950-3100\t20\t11\t0.95\tNA\t13\t70\tViruses;Caudoviricetes\n")
        table = run(root)
        det = [r for r in table if r["called"] == "detected_target_overlap"]
        assert det and det[0]["coordinates"] == "950-3100"
        assert det[0]["overlap_with_element"].startswith("2001 bp")
        print("A (geNomad target detected)        OK")
    finally:
        shutil.rmtree(root, ignore_errors=True)

    root = fresh()
    try:
        write(os.path.join(root, "genomad", "smbu06", "run_virus_summary.tsv"),
              GSM_HDR + "Chr|provirus_50000_52000\t2000\tProvirus\t50000-52000\t10\t11\t0.9\tNA\t5\t20\tV;G\n")
        write(os.path.join(root, "genomad", "element", "run_virus_summary.tsv"),
              GSM_HDR + "element\t2000\tNo terminal repeats\tNA\t10\t11\t0.9\tNA\t5\t20\tV;G\n")
        table = run(root)
        assert [r for r in table if r["called"] == "no_target_overlap"]
        assert [r for r in table if r["called"].startswith("virus (sequence-level")]
        print("D (geNomad no target overlap)      OK")
    finally:
        shutil.rmtree(root, ignore_errors=True)

    root = fresh()
    try:
        write(os.path.join(root, "genomad", "smbu06", "run_virus_summary.tsv"), GSM_HDR)  # header only
        table = run(root)
        empt = [r for r in table if r["called"] == "run_failed or empty output"]
        assert empt, f"empty summary not reported: {table}"
        print("E (geNomad empty summary)          OK")
    finally:
        shutil.rmtree(root, ignore_errors=True)

    # ---------------- PhiSpy coordinate parsing ----------------
    hits = phispy_case([], ["detected_target_overlap"],
                       "pp1\tChromosome1\t1200\t2800\n")   # headerless, over target
    assert hits[0]["coordinates"] == "1200-2800", hits[0]
    print("F (headerless single, target)      OK")

    hits = phispy_case([], ["retained call elsewhere"],
                       "pp1\tPlasmid1\t110154\t128787\n", cands=False)
    assert any("110154-128787" in h["coordinates"] for h in hits), hits
    print("G (headerless single, elsewhere)   OK")

    hits = phispy_case([], ["detected_target_overlap", "retained call elsewhere"],
                       "pp1\tChromosome1\t1200\t2800\npp2\tPlasmid1\t100\t200\n")
    assert len(hits) >= 2
    print("H (multiple calls)                 OK")

    hits = phispy_case([], ["detected_target_overlap"],
                       "pp\tcontig\tstart\tstop\tatts\npp1\tChromosome1\t1200\t2800\t\n")  # WITH header
    assert hits[0]["coordinates"] == "1200-2800", hits[0]
    print("I (header tolerated)               OK")

    # B/J: the element-region candidate is present in the log but the (empty,
    # completed) coordinates file retains nothing -> reported as "candidate dropped"
    hits = phispy_case([], ["candidate dropped"], "", cands=True)
    assert "No genes were identified" in hits[0]["coordinates"], hits[0]
    print("B/J (candidate dropped, completed) OK")

    hits = phispy_case([], ["unparseable"], "pp1\tChromosome1\tBROKEN\t2800\n")
    assert "line 1" in hits[0]["notes"], hits[0]
    print("K (malformed line reported)        OK")

    # ---------------- PhiSpy run-state gates ----------------
    hits = phispy_case([], ["incomplete"], "pp1\tPlasmid1\t100\t200\n", done=False)
    assert not any("no_target_overlap" in h["called"] for h in hits), hits  # NOT a success-negative
    print("L (truncated log)                  OK")

    hits = phispy_case([], ["no candidates logged"], "", no_cand_table=True)
    print("M (completed, no candidates)       OK")

    hits = phispy_case([], ["coordinates file missing"], None, cands=True, write_coords=False)
    notes = hits[0]["notes"]
    assert "cannot be inferred" in notes, notes  # retention NOT inferred as dropped
    print("N (candidates but coords missing)  OK")

    root = fresh()
    try:
        write(os.path.join(root, "phispy", "phispy_smbu06.log"),
              "Traceback (most recent call last):\nModuleNotFoundError: x\n")
        table = run(root)
        rf = [r for r in table if r["tool"].startswith("PhiSpy") and r["called"] == "run_failed"]
        assert rf, table
        print("O (explicit failure)               OK")
    finally:
        shutil.rmtree(root, ignore_errors=True)

    hits = phispy_case([], ["incomplete"], "pp1\tChromosome1\t1200\t2800\n", done=False)
    assert "inconsistent" in hits[0]["notes"] and "cannot be confirmed" in hits[0]["notes"], hits[0]
    print("P (candidates+rows w/o completion) OK")

    print("ALL SCENARIOS PASSED")


if __name__ == "__main__":
    main()
