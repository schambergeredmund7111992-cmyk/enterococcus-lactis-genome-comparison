#!/usr/bin/env python3
"""Run geNomad end-to-end on the project inputs (independent element re-annotation).

Safety model for the output directory (OUTPUT):
  * OUTPUT does not exist            -> run normally.
  * OUTPUT exists, no flag           -> refuse; suggest --restart (default-safe,
                                        geNomad reuses/overwrites intermediate files
                                        without deleting the directory) or an explicit
                                        --overwrite.
  * OUTPUT exists, --restart         -> pass --restart to geNomad (nothing is deleted
                                        by this script).
  * OUTPUT exists, --overwrite       -> delete ONLY after all safety checks pass:
        - OUTPUT (fully resolved) must not be a drive root, the user home, the
          repository root or any ancestor of it;
        - OUTPUT must not equal, contain, or live inside the input file, its
          directory, or the database directory (checked on resolved, real paths);
        - the path as given must contain no symbolic link or Windows junction
          component (otherwise the deletion target is ambiguous -> refuse);
        - the existing directory must look like a geNomad output of this tool
          (*_annotate.json / *_summary/ / *_annotate.log present), or be empty;
          anything else is not deleted.

Usage:
    python run_genomad.py --input <fasta> --out <out_dir> --db <genomad_db> \
        [--restart | --overwrite] [--splits 4] [--threads 6]

Notes
-----
* geNomad v1.12.0 with database v1.9 (2025-02-18) and MMseqs2 18.8cc5c were used
  for the published results (all four inputs completed). A previous bioconda
  pairing (geNomad 1.7.0 + MMseqs2 14.7e284) failed at the marker-database
  prefilter step; see docs/element_tool_validation_methods.md.
* Assembly/element sequences are provider deliveries and are not distributed in
  this repository (see input_manifest.tsv).
"""
import argparse
import os
import shutil
import subprocess
import sys


def resolve(p):
    """Fully resolved absolute path (resolves symlinks and Windows junctions)."""
    return os.path.realpath(os.path.abspath(p))


def path_has_link_component(p):
    """True if any component of the path as given is a symlink or junction."""
    p = os.path.abspath(p)
    parts = []
    while True:
        head, tail = os.path.split(p)
        if tail:
            parts.append(p)
            p = head
        else:
            break
    for q in parts:
        if os.path.islink(q) or (hasattr(os.path, "isjunction") and os.path.isjunction(q)):
            return True
    return False


def find_repo_root(start):
    p = resolve(start)
    while True:
        if os.path.isdir(os.path.join(p, ".git")):
            return p
        parent = os.path.dirname(p)
        if parent == p:
            return None
        p = parent


def is_ancestor_or_equal(a, b):
    """True if resolved path a is an ancestor of (or equal to) resolved path b."""
    a, b = resolve(a), resolve(b)
    return b == a or b.startswith(a + os.sep)


def looks_like_genomad_output(d):
    import glob
    if glob.glob(os.path.join(d, "*_annotate.json")) or glob.glob(os.path.join(d, "*_annotate.log")):
        return True
    for sub in glob.glob(os.path.join(d, "*")):
        if os.path.isdir(sub) and (sub.endswith("_summary") or sub.endswith("_annotate")):
            return True
    return False


def refuse(msg):
    sys.exit(f"REFUSED: {msg}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="input FASTA (element, assembly or component)")
    ap.add_argument("--out", required=True, help="output directory (geNomad sample dir)")
    ap.add_argument("--db", required=True, help="geNomad database directory")
    ap.add_argument("--splits", type=int, default=4)
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--genomad", default="genomad", help="genomad executable (default: on PATH)")
    ap.add_argument("--mmseqs-bin", default=None, help="directory containing mmseqs (optional)")
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--restart", action="store_true",
                      help="reuse/overwrite intermediate files in an existing OUTPUT (nothing deleted by this script)")
    mode.add_argument("--overwrite", action="store_true",
                      help="explicitly delete an existing OUTPUT directory after safety checks")
    args = ap.parse_args()

    if not os.path.exists(args.input):
        sys.exit(f"input not found: {args.input}")
    if not os.path.isdir(args.db):
        sys.exit(f"database dir not found: {args.db}")
    if args.mmseqs_bin:
        os.environ["PATH"] = args.mmseqs_bin + os.pathsep + os.environ["PATH"]

    out = args.out
    extra = []
    if os.path.exists(out):
        if not os.path.isdir(out):
            refuse(f"OUTPUT exists and is not a directory: {out}")
        if args.overwrite:
            r_out = resolve(out)
            # --- safety checks before any deletion ---
            if os.path.dirname(r_out) == r_out:
                refuse(f"OUTPUT resolves to a drive root: {r_out}")
            home = resolve(os.path.expanduser("~"))
            if r_out == home:
                refuse(f"OUTPUT resolves to the user home: {r_out}")
            repo = find_repo_root(os.path.dirname(os.path.abspath(__file__)))
            if repo and is_ancestor_or_equal(r_out, repo):
                refuse(f"OUTPUT is the repository root or an ancestor of it: {r_out}")
            if path_has_link_component(out):
                refuse("the OUTPUT path contains a symlink/junction component; "
                       "pass the real directory path explicitly")
            for protected, what in ((os.path.dirname(resolve(args.input)), "the input directory"),
                                    (args.db, "the database directory")):
                r_p = resolve(protected)
                if r_out == r_p:
                    refuse(f"OUTPUT is {what}: {r_p}")
                if is_ancestor_or_equal(r_out, r_p):
                    refuse(f"OUTPUT would delete {what} (it is an ancestor of it): {r_out} ⊇ {r_p}")
                if is_ancestor_or_equal(r_p, r_out):
                    refuse(f"OUTPUT lies inside {what}: {r_out}")
            contents = os.listdir(out)
            if contents and not looks_like_genomad_output(out):
                refuse("OUTPUT is not empty and does not look like a geNomad output of this tool; "
                       "refusing to delete it")
            print(f"removing existing geNomad output directory (explicit --overwrite): {out}")
            shutil.rmtree(out)
        elif args.restart:
            print(f"reusing existing output directory with --restart (nothing deleted): {out}")
            extra = ["--restart"]
        else:
            refuse(f"OUTPUT exists: {out}\n"
                   "  Use --restart to reuse/overwrite intermediate files without deleting, "
                   "or --overwrite to explicitly delete and re-run.")
    elif args.restart:
        print("note: --restart given but OUTPUT does not exist; running normally")

    cmd = [args.genomad, "end-to-end", "--cleanup", "--splits", str(args.splits),
           "--threads", str(args.threads)] + extra + [args.input, out, args.db]
    print("RUN:", " ".join(cmd), flush=True)
    r = subprocess.run(cmd)
    sys.exit(r.returncode)


if __name__ == "__main__":
    main()
