#!/usr/bin/env python3
"""Run PhiSpy on a GenBank record (second, methodologically independent tool).

Usage:
    python run_phispy.py --genbank <assembly.gb> --out <out_dir> [--threads 4] \
        [--phispy <phispy executable>]

Notes
-----
* PhiSpy v5.0.10 (PyPI package; Akhter et al. 2012) was run with default
  settings and WITHOUT the optional --phmms module (hmmsearch not required for
  the default pipeline).
* Input is a GenBank file with provider annotation. Used for the published run:
  the complete smbu06 and smbu08 assemblies (GenBank). Sequence files are not
  distributed in this repository.
* Candidate regions are ranked by the classifier and then filtered: in the
  v5.0.10 evaluation code a candidate is kept only when at least one of its
  genes is classified as a phage gene (default threshold 1). Both project
  genomes produced no retained call for the element region; see
  results/tables/element_tool_validation/phispy/ and
  docs/element_tool_validation_methods.md.
"""
import argparse
import os
import subprocess
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genbank", required=True, help="input GenBank file")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--phispy", default="phispy", help="phispy executable (default: on PATH)")
    args = ap.parse_args()

    if not os.path.exists(args.genbank):
        sys.exit(f"genbank not found: {args.genbank}")
    os.makedirs(args.out, exist_ok=True)
    cmd = [args.phispy, args.genbank, "-o", args.out, "--threads", str(args.threads)]
    print("RUN:", " ".join(cmd), flush=True)
    r = subprocess.run(cmd)
    sys.exit(r.returncode)


if __name__ == "__main__":
    main()
