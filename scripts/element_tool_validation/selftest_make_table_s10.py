#!/usr/bin/env python3
"""Self-test for make_table_s10.py using small synthetic fixtures.

Covers four scenarios:
  A. target detected      - geNomad provirus interval overlapping the target
  B. candidate dropped    - PhiSpy ranked candidate overlapping the target, not retained
  C. input missing        - an expected run directory absent
  D. no target overlap    - calls present but disjoint from the target interval
Each scenario is generated in a temporary directory; make_table_s10.py is run as a
subprocess and its output table is checked. No real study data are touched.

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


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="\n") as f:
        f.write(text)


def run(root, tag_target=("Chromosome1", 1000, 3000)):
    out = os.path.join(root, "s10.tsv")
    cmd = [sys.executable, SCRIPT,
           "--genomad-root", os.path.join(root, "genomad"),
           "--phispy-root", os.path.join(root, "phispy"),
           "--target-chromosome", tag_target[0],
           "--target-start", str(tag_target[1]),
           "--target-end", str(tag_target[2]),
           "--out", out]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("STDERR:", r.stderr)
        raise SystemExit("script failed")
    with open(out, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def scenario_A(root):
    """geNomad provirus overlapping the target -> detected_target_overlap."""
    write(os.path.join(root, "genomad", "smbu06", "run_virus_summary.tsv"),
          "seq_name\tlength\ttopology\tcoordinates\tn_genes\tgenetic_code\tvirus_score\tfdr\tn_hallmarks\tmarker_enrichment\ttaxonomy\n"
          "Chr|provirus_950_3100\t2150\tProvirus\t950-3100\t20\t11\t0.95\tNA\t13\t70\tViruses;Caudoviricetes\n")
    write(os.path.join(root, "genomad", "smbu06", "run_marker_classification.tsv"),
          "seq_name\tchromosome_score\tplasmid_score\tvirus_score\nChr\t0.1\t0.1\t0.9\n")
    write(os.path.join(root, "genomad", "smbu06", "run_nn_classification.tsv"),
          "seq_name\tchromosome_score\tplasmid_score\tvirus_score\nChr\t0.2\t0.2\t0.6\n")
    rows = run(root)
    det = [r for r in rows if r["called"] == "detected_target_overlap"]
    assert det and det[0]["coordinates"] == "950-3100", f"A failed: {rows}"
    assert det[0]["overlap_with_element"].startswith("2001 bp"), f"A overlap wrong: {det[0]}"
    print("A (target detected)               OK")


def scenario_B(root):
    """PhiSpy candidate overlapping target, dropped -> no retained call; message parsed."""
    write(os.path.join(root, "phispy", "smbu06.log"),
          "Potential prophages (sorted highest to lowest)\n"
          "Contig\tStart\tStop\tNumber of potential genes\tStatus\n"
          "Chromosome1\t1200\t2800\t24\tDropped. No genes were identified as phage genes\n"
          "Plasmid1\t200\t400\t2\tDropped. Region too small (Not enough genes)\n"
          "Done!!!\n")
    write(os.path.join(root, "phispy", "prophage_coordinates_smbu06.tsv"),
          "pp\tcontig\tstart\tstop\n")
    rows = run(root)
    hit = [r for r in rows if r["tool"].startswith("PhiSpy") and "smbu06" in r["input"]]
    assert hit, f"B failed: {rows}"
    assert "candidate dropped" in hit[0]["called"], f"B status wrong: {hit[0]}"
    assert "1264149" not in hit[0]["coordinates"] and "1200-2800" in hit[0]["coordinates"], f"B coords: {hit[0]}"
    assert "No genes were identified" in hit[0]["coordinates"], f"B message: {hit[0]}"
    print("B (candidate dropped)             OK")


def scenario_C(root):
    """Absent geNomad run dir -> input_missing; no fabricated result rows."""
    rows = run(root)
    miss = [r for r in rows if r["called"] == "input_missing"]
    assert miss, f"C failed: {rows}"
    assert all(r["coordinates"] == "" and r["score"] == "" for r in miss), "C fabricated values"
    print("C (input missing)                 OK")


def scenario_D(root):
    """Calls disjoint from target -> no_target_overlap; sequence-level input flagged."""
    write(os.path.join(root, "genomad", "smbu06", "run_virus_summary.tsv"),
          "seq_name\tlength\ttopology\tcoordinates\tn_genes\tgenetic_code\tvirus_score\tfdr\tn_hallmarks\tmarker_enrichment\ttaxonomy\n"
          "Chr|provirus_50000_52000\t2000\tProvirus\t50000-52000\t10\t11\t0.9\tNA\t5\t20\tViruses;Caudoviricetes\n")
    write(os.path.join(root, "genomad", "element", "run_virus_summary.tsv"),
          "seq_name\tlength\ttopology\tcoordinates\tn_genes\tgenetic_code\tvirus_score\tfdr\tn_hallmarks\tmarker_enrichment\ttaxonomy\n"
          "element\t2000\tNo terminal repeats\tNA\t10\t11\t0.9\tNA\t5\t20\tViruses;Caudoviricetes\n")
    rows = run(root)
    nto = [r for r in rows if r["called"] == "no_target_overlap"]
    seq = [r for r in rows if r["called"].startswith("virus (sequence-level")]
    assert nto, f"D failed: {rows}"
    assert seq, f"D seq-level failed: {rows}"
    assert nto[0]["overlap_with_element"] == "none"
    print("D (no target overlap)             OK")


def main():
    print("self-test: make_table_s10.py")
    # B and C share fixtures: run C first (missing dirs), then add B fixtures and re-run
    root = tempfile.mkdtemp(prefix="s10test_")
    try:
        scenario_C(root)                 # nothing exists yet -> input_missing
        scenario_B(root)                 # add phispy
        shutil.rmtree(os.path.join(root, "genomad"), ignore_errors=True)
        scenario_A(root)                 # add geNomad with overlapping provirus
        shutil.rmtree(os.path.join(root, "genomad"), ignore_errors=True)
        scenario_D(root)                 # disjoint + sequence-level
        print("ALL SCENARIOS PASSED")
    finally:
        shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    main()
