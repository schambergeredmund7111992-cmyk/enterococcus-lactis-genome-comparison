#!/usr/bin/env python3
"""Manuscript Figure 5 (V9): candidate bacteriocin-associated loci, structural detail.

Usage:
    python fig5_loci.py --out-prefix <path/Fig5_manuscript>

Replaces the earlier Figure 5, whose panel (a) duplicated Figure 4's backbone
view. Panels show the two candidate loci with provider-annotation coordinates,
the 10-bp spacing to the candidate immunity factor, and the read-supported
upstream in-frame stop codons (no complete precursor ORF; Section 3.5).

Coordinates (shared plasmid, 1-based, provider annotation):
    GL002641 enterocin P-like        36,073-36,237  (+)
    GL002642 candidate immunity      36,247-36,513  (+), 10 bp downstream
    GL002679 hiracin JM79-like       69,485-69,688  (+), truncated homologue
Upstream in-frame stops: TGA 19 codons upstream of GL002641; TAG 3 codons
upstream of GL002679 (read-supported: 97.4% of 424 and 97.8% of 499 Nanopore
reads). No sequence data are required or distributed.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ORANGE = "#D95F02"
GREEN = "#3A9E6E"
RED = "#C0392B"
GREY = "#8C8C8C"
TNR = "Times New Roman"
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = [TNR, "DejaVu Serif"]


def gene(ax, s, e, col, label):
    w = e - s
    head = min(w * 0.25, 25)
    ax.add_patch(plt.Polygon([(s, 0.35), (e - head, 0.35), (e, 0.5), (e - head, 0.65), (s, 0.65)],
                             closed=True, facecolor=col, edgecolor="#333333", lw=0.4, zorder=3))
    ax.text((s + e) / 2, 0.72, label, ha="center", va="bottom", fontsize=5.6, color=col, zorder=3)


def stopmark(ax, pos, txt):
    ax.plot([pos, pos], [0.30, 0.70], color=RED, lw=1.1, zorder=4)
    ax.text(pos, 0.80, txt, ha="center", va="bottom", fontsize=5.3, color=RED, zorder=4)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    fig = plt.figure(figsize=(7.62, 4.03))
    gs = fig.add_gridspec(2, 1, hspace=0.62, left=0.05, right=0.97, top=0.90, bottom=0.09)
    fig.suptitle("Candidate bacteriocin-associated loci on the shared plasmid (structural detail)",
                 fontsize=8.4, fontweight="bold", x=0.05, ha="left", y=0.97)

    ax = fig.add_subplot(gs[0])
    ax.set_xlim(35940, 36760); ax.set_ylim(0.05, 1.30); ax.axis("off")
    ax.text(-0.01, 1.00, "a", transform=ax.transAxes, fontsize=10, fontweight="bold")
    ax.text(0.05, 1.04, "enterocin P-like locus (shared plasmid 35,940-36,760)",
            transform=ax.transAxes, fontsize=6.6, fontweight="bold")
    ax.text(0.055, 0.90, "no complete precursor ORF is encoded\n(read-supported in-frame stop)",
            transform=ax.transAxes, ha="left", va="top", fontsize=5.4, color="#444444")
    ax.plot([35960, 36740], [0.5, 0.5], color=GREY, lw=0.8, zorder=1)
    gene(ax, 36073, 36237, ORANGE, "enterocin P-like\n(mature-peptide region)")
    gene(ax, 36247, 36513, GREEN, "candidate immunity factor")
    stopmark(ax, 36016, "in-frame TGA stop\n19 codons upstream")
    ax.annotate("", xy=(36257, 0.22), xytext=(36227, 0.22),
                arrowprops=dict(arrowstyle="<->", lw=0.6, color="#333333"))
    ax.text(36242, 0.10, "10 bp", ha="center", fontsize=5.2)
    ax.text(0.985, 0.42, "...", transform=ax.transAxes, fontsize=7, color=GREY)

    ax = fig.add_subplot(gs[1])
    ax.set_xlim(69260, 70060); ax.set_ylim(0.05, 1.30); ax.axis("off")
    ax.text(-0.01, 1.00, "b", transform=ax.transAxes, fontsize=10, fontweight="bold")
    ax.text(0.05, 1.04, "hiracin JM79-like locus (shared plasmid 69,260-70,060)",
            transform=ax.transAxes, fontsize=6.6, fontweight="bold")
    ax.text(0.055, 0.90, "truncated relative to its reference precursor;\nno complete precursor ORF",
            transform=ax.transAxes, ha="left", va="top", fontsize=5.4, color="#444444")
    ax.plot([69280, 70040], [0.5, 0.5], color=GREY, lw=0.8, zorder=1)
    gene(ax, 69485, 69688, RED, "hiracin JM79-like\n(truncated homologue)")
    stopmark(ax, 69476, "in-frame TAG stop\n3 codons upstream")
    ax.text(0.965, 0.42, "...", transform=ax.transAxes, fontsize=7, color=GREY)

    fig.text(0.05, 0.012,
             "Coordinates: provider annotation of the shared plasmid. Identities and coverage: Section 3.5; "
             "stop-codon support: read pileups (97.4% of 424 and 97.8% of 499 Nanopore reads).",
             fontsize=5.0, color="#444444")

    os.makedirs(os.path.dirname(os.path.abspath(args.out_prefix)), exist_ok=True)
    for ext in (".png", ".pdf", ".svg"):
        fig.savefig(args.out_prefix + ext, dpi=600, bbox_inches="tight")
    print("saved", args.out_prefix + ".{png,pdf,svg}")


if __name__ == "__main__":
    main()
