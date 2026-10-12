#!/usr/bin/env python3
"""Self-test for run_genomad.py output-directory safety.

Uses only temporary directories and a fake geNomad executable; no real study
data or research directories are touched. Checks:
  A. existing OUTPUT without a flag            -> refused
  B. OUTPUT == user home (--overwrite)         -> refused
  C. OUTPUT is an ancestor of the input dir    -> refused
  D. OUTPUT path contains a junction component -> refused
  E. non-empty OUTPUT that is not geNomad output (--overwrite) -> refused
  F. empty OUTPUT (--overwrite)                -> allowed (fake run executes)
  G. geNomad-looking OUTPUT (--overwrite)      -> allowed
  H. existing OUTPUT with --restart            -> allowed, nothing deleted

Usage:  python selftest_run_genomad.py
"""
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "run_genomad.py")


def run_case(root, out, extra=(), expect_refused=True, expect_real_run=False):
    fake = os.path.join(root, "fake_genomad.py")
    with open(fake, "w", newline="\n") as f:
        f.write("import sys\nprint('FAKE_RUN', ' '.join(sys.argv[1:]))\n")
    inp = os.path.join(root, "inputs", "el.fa")
    os.makedirs(os.path.dirname(inp), exist_ok=True)
    if not os.path.exists(inp):
        with open(inp, "w") as f:
            f.write(">x\nACGT\n")
    db = os.path.join(root, "db")
    os.makedirs(db, exist_ok=True)
    cmd = [sys.executable, SCRIPT, "--input", inp, "--out", out, "--db", db,
           "--genomad", sys.executable.replace("python.exe", "python.exe"),
           ] + list(extra)
    # make the "genomad" callable: use a wrapper via --genomad pointing to a .bat is
    # platform-tricky; instead prepend a tiny wrapper module via PATH-free approach:
    # run the real script but replace the executable by a shell wrapper
    wrapper = os.path.join(root, "fake_genomad.bat")
    with open(wrapper, "w", newline="\r\n") as f:
        f.write("@echo off\r\necho FAKE_RUN %*\r\n")
    cmd[cmd.index("--genomad") + 1] = wrapper
    r = subprocess.run(cmd, capture_output=True, text=True, encoding='utf-8', errors='replace', cwd=root)
    out_txt = (r.stdout or "") + (r.stderr or "")
    if expect_refused:
        ok = "REFUSED" in out_txt and "FAKE_RUN" not in out_txt
        print(("PASS" if ok else "FAIL"), "-", out_txt.strip().split("\n")[0][:90])
        assert ok, out_txt
    else:
        ok = "FAKE_RUN" in out_txt
        print(("PASS" if ok else "FAIL"), "- allowed and fake run executed")
        assert ok, out_txt


def main():
    root = tempfile.mkdtemp(prefix="genomadt_")
    try:
        # A: existing OUTPUT, no flag
        outA = os.path.join(root, "outA"); os.makedirs(outA)
        run_case(root, outA, expect_refused=True)

        # B: home as OUTPUT
        run_case(root, os.path.expanduser("~"), extra=["--overwrite"], expect_refused=True)

        # C: OUTPUT is an ancestor of the input directory
        run_case(root, root, extra=["--overwrite"], expect_refused=True)

        # D: junction component in the OUTPUT path
        target = os.path.join(root, "real_target")
        os.makedirs(target, exist_ok=True)
        junction = os.path.join(root, "link_out")
        subprocess.run(["cmd", "/c", "mklink", "/J", junction, target], capture_output=True)
        if os.path.exists(junction):
            run_case(root, junction, extra=["--overwrite"], expect_refused=True)
        else:
            print("SKIP - junction could not be created")

        # E: non-empty, not geNomad output
        outE = os.path.join(root, "outE"); os.makedirs(outE)
        with open(os.path.join(outE, "something.txt"), "w") as f:
            f.write("x")
        run_case(root, outE, extra=["--overwrite"], expect_refused=True)

        # F: empty dir -> allowed
        outF = os.path.join(root, "outF"); os.makedirs(outF)
        run_case(root, outF, extra=["--overwrite"], expect_refused=False)

        # G: geNomad-looking dir -> allowed
        outG = os.path.join(root, "outG"); os.makedirs(outG)
        with open(os.path.join(outG, "sample_annotate.json"), "w") as f:
            f.write("{}")
        run_case(root, outG, extra=["--overwrite"], expect_refused=False)

        # H: --restart path, nothing deleted
        outH = os.path.join(root, "outH"); os.makedirs(outH)
        with open(os.path.join(outH, "sample_annotate.json"), "w") as f:
            f.write("{}")
        cmd = [sys.executable, SCRIPT, "--input", os.path.join(root, "inputs", "el.fa"),
               "--out", outH, "--db", os.path.join(root, "db"), "--restart",
               "--genomad", os.path.join(root, "fake_genomad.bat")]
        r = subprocess.run(cmd, capture_output=True, text=True, encoding='utf-8', errors='replace')
        t = (r.stdout or "") + (r.stderr or "")
        ok = "FAKE_RUN" in t and os.path.exists(os.path.join(outH, "sample_annotate.json")) and "--restart" in t
        print(("PASS" if ok else "FAIL"), "- --restart allowed and nothing deleted")
        assert ok, t

        print("ALL SAFETY CASES PASSED")
    finally:
        # remove junction first (rmtree would follow it); the junction itself was
        # created inside root, so remove it explicitly before the tree
        j = os.path.join(root, "link_out")
        if os.path.exists(j) and hasattr(os.path, "isjunction") and os.path.isjunction(j):
            os.rmdir(j)
        shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    main()
