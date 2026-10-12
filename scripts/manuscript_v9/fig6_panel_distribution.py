#!/usr/bin/env python3
"""Manuscript Figure 6 (V9): public-panel backbone distribution and candidate-gene
co-occurrence. All counts are computed from the input tables; none are hardcoded.

Usage:
    python fig6_panel_distribution.py \
        --distribution ../results/tables/supplementary/TableS4_homology_thresholds/pl1_distribution.tsv \
        --presence     ../results/tables/supplementary/TableS4_homology_thresholds/gene_presence_matrix.tsv \
        --panel-metadata <optional: panel_metadata.tsv, used only to name the six
                           panel genomes without detectable backbone alignments> \
        --out-prefix <path/Fig6_manuscript>

Inclusion rules and denominators (from the pipeline that produced the inputs):
* panel = 159 deduplicated public genomes (accessions in Table S3 of the
  repository / manuscript Table S4); all are complete assemblies.
* (a) backbone distribution: a genome is classified only if it has >=1 detectable
  BLASTN alignment to the 128,837-bp plasmid (153 of 159); classes: near-identical
  (>=95% colinear-chain coverage and >=99% length-weighted identity), partial
  (>=30% coverage), fragment (>=5% coverage), trace (<5% coverage). The six
  genomes without any alignment (all Lactococcus controls) are excluded from the
  panel and their identities are printed by this script.
* (b) candidate-gene co-occurrence: computed over all 159 panel genomes from the
  presence matrix at the 70%/70% threshold (>=70% identity, >=70% query
  coverage); the Fisher exact test in the manuscript uses the same sample set.
"""
import argparse
import csv
import os
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TNR = "Times New Roman"
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = [TNR, "DejaVu Serif"]

CLASS_ORDER = ["full_length_near_identical", "partial_backbone", "fragment", "trace"]
CLASS_LABELS = {"full_length_near_identical": "near-identical\n(>95% cov, >99% id)",
                "partial_backbone": "partial\n(>30% cov)",
                "fragment": "fragment\n(>5% cov)",
                "trace": "trace\n(<5% cov)"}
CLASS_COLORS = {"full_length_near_identical": "#2B6CB0", "partial_backbone": "#7BA7D7",
                "fragment": "#B9CDE5", "trace": "#D9D9D9"}
# NOTE: combination keys are sorted tuples (e.g. ("H","I","P")) to match the
# computation below; labels are written in the conventional P/I/H order.
GENE_ORDER = [("H", "I", "P"), ("I", "P"), ("P",), ("H",)]
GENE_LABELS = {("H", "I", "P"): "P + I + H", ("I", "P"): "P + I",
               ("P",): "P only", ("H",): "H only"}
GENE_COLORS = {("H", "I", "P"): "#C0392B", ("I", "P"): "#E67E22",
               ("P",): "#C0392B", ("H",): "#BDBDBD"}


def read_tsv(p):
    with open(p, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--distribution", required=True)
    ap.add_argument("--presence", required=True)
    ap.add_argument("--panel-metadata", default=None)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    # ---- (a) backbone distribution ----
    dist_rows = read_tsv(args.distribution)
    class_counts = Counter(r["class"] for r in dist_rows)
    with_alignment = {r["accession"] for r in dist_rows}

    # ---- (b) co-occurrence over the full panel ----
    pres_rows = read_tsv(args.presence)
    per_genome = defaultdict(dict)
    for r in pres_rows:
        per_genome[r["accession"]][r["gene"]] = r["presence_70_70"]
    name_map = {"enterocinP_like": "P", "sakacinA_immunity": "I", "hiracinJM79_like": "H"}
    combo_counts = Counter()
    for acc, gd in per_genome.items():
        combo = tuple(sorted(name_map[g] for g, v in gd.items() if v == "yes"))
        combo_counts[combo] += 1

    panel_n = len(per_genome)
    no_align = sorted(set(per_genome) - with_alignment)
    print(f"panel genomes: {panel_n}; with backbone alignment: {len(with_alignment)}")
    print("class counts:", {k: class_counts.get(k, 0) for k in CLASS_ORDER})
    print("no-alignment genomes:", no_align)
    if args.panel_metadata and os.path.exists(args.panel_metadata):
        meta = {r["accession"]: r for r in read_tsv(args.panel_metadata)}
        for a in no_align:
            m = meta.get(a, {})
            print(f"   {a}: {m.get('organism', '?')} | {m.get('strain', '?')}")

    # ---- figure ----
    fig = plt.figure(figsize=(9.5, 5.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0], wspace=0.28,
                          left=0.075, right=0.98, top=0.88, bottom=0.12)

    ax = fig.add_subplot(gs[0])
    vals = [class_counts.get(k, 0) for k in CLASS_ORDER]
    bars = ax.bar(range(len(CLASS_ORDER)), vals,
                  color=[CLASS_COLORS[k] for k in CLASS_ORDER], edgecolor="#333333", lw=0.5)
    for x, v in zip(range(len(CLASS_ORDER)), vals):
        ax.text(x, v + max(vals) * 0.02, str(v), ha="center", fontsize=8)
    ax.set_xticks(range(len(CLASS_ORDER)))
    ax.set_xticklabels([CLASS_LABELS[k] for k in CLASS_ORDER], fontsize=7)
    ax.set_ylabel("genomes (n)", fontsize=8)
    ax.set_title("a  Plasmid-backbone distribution", fontsize=9.5, fontweight="bold", loc="left")
    ax.text(0.98, 0.95, f"n = {len(with_alignment)}", transform=ax.transAxes,
            ha="right", va="top", fontsize=9)
    ax.tick_params(axis="y", labelsize=7)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    ax2 = fig.add_subplot(gs[1])
    cvals = [combo_counts.get(k, 0) for k in GENE_ORDER]
    ax2.barh(range(len(GENE_ORDER))[::-1], cvals,
             color=[GENE_COLORS[k] for k in GENE_ORDER], edgecolor="#333333", lw=0.5)
    for y, v in zip(range(len(GENE_ORDER))[::-1], cvals):
        ax2.text(v + max(cvals) * 0.02, y, str(v), va="center", fontsize=8)
    ax2.set_yticks(range(len(GENE_ORDER))[::-1])
    ax2.set_yticklabels([GENE_LABELS[k] for k in GENE_ORDER], fontsize=7.5)
    ax2.set_xlabel("genomes (n)", fontsize=8)
    ax2.set_title("b  Candidate-gene co-occurrence (70%/70%)", fontsize=9.5,
                  fontweight="bold", loc="left")
    ax2.tick_params(axis="x", labelsize=7)
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)

    os.makedirs(os.path.dirname(os.path.abspath(args.out_prefix)), exist_ok=True)
    for ext in (".png", ".pdf", ".svg"):
        fig.savefig(args.out_prefix + ext, dpi=600, bbox_inches="tight")
    print("saved", args.out_prefix + ".{png,pdf,svg}")


if __name__ == "__main__":
    main()
