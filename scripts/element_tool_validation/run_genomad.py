#!/usr/bin/env python3
"""Run geNomad end-to-end on the project inputs (independent element re-annotation).

Usage:
    python run_genomad.py --input <assembly_or_element.fasta> --out <out_dir> \
        --db <genomad_db_dir> [--splits 4] [--threads 6]

Notes
-----
* geNomad v1.12.0 with database v1.9 (2025-02-18 release) and MMseqs2 18.8cc5c
  were used for the published results. Install: conda/mamba install -c bioconda
  genomad; download the database with `genomad download-database <dir>`.
* The published runs used the element sequence (attL-attR, 38,865 bp), each
  complete assembly (both replicons) and the 77,437-bp circular component as
  separate inputs. Assembly sequences are not distributed in this repository;
  provide your own copy of the deliveries (see input_manifest.tsv).
* A previous bioconda pairing (geNomad 1.7.0 with MMseqs2 14.7e284) failed at
  the IndexTable memory-allocation step for the marker database; geNomad 1.12.0
  + MMseqs2 18.8cc5c completed all four inputs (see docs/element_tool_validation_methods.md).
"""
import argparse
import os
import shutil
import subprocess
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="input FASTA (element, assembly or component)")
    ap.add_argument("--out", required=True, help="output directory (geNomad sample dir)")
    ap.add_argument("--db", required=True, help="geNomad database directory")
    ap.add_argument("--splits", type=int, default=4)
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--genomad", default="genomad", help="genomad executable (default: on PATH)")
    ap.add_argument("--mmseqs-bin", default=None, help="directory containing mmseqs (optional)")
    args = ap.parse_args()

    if not os.path.exists(args.input):
        sys.exit(f"input not found: {args.input}")
    if not os.path.isdir(args.db):
        sys.exit(f"database dir not found: {args.db}")
    if args.mmseqs_bin:
        os.environ["PATH"] = args.mmseqs_bin + os.pathsep + os.environ["PATH"]

    out = args.out
    if os.path.exists(out):
        shutil.rmtree(out)
    cmd = [args.genomad, "end-to-end", "--cleanup", "--splits", str(args.splits),
           "--threads", str(args.threads), args.input, out, args.db]
    print("RUN:", " ".join(cmd), flush=True)
    r = subprocess.run(cmd)
    sys.exit(r.returncode)


if __name__ == "__main__":
    main()
