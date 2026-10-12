"""Manuscript Supplementary Figure S6 (V9): synteny comparison.

Usage:
    python plot_synteny.py --data-dir <dir with homology_links.tsv, all_region_genes.tsv,
        element_gene_meta.tsv> --out-prefix <path/FigS6_synteny>

The input tables are published in results/tables/element_tool_validation/synteny/.

Key changes vs v1:
 - Every row is trimmed to its element-homologous span (+-2 kb) so rows align;
   PBSX shows its complete 25,226-bp region; skin its complete 50-kb region.
 - Element row coloured from its own 51-CDS module assignment.
 - PBSX / skin rows coloured from their own GenBank /product keyword assignment,
   demonstrating module architecture similarity WITHOUT implying sequence homology.
 - Links only from real TBLASTN hits (pident>=30, qcov>=50%, gene overlap>=50%).
"""
import csv
import os
import re
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import matplotlib.patches as mpatches

import argparse as _argparse
_ap = _argparse.ArgumentParser()
_ap.add_argument("--data-dir", required=True, help="dir with homology_links.tsv, all_region_genes.tsv, element_gene_meta.tsv")
_ap.add_argument("--out-prefix", required=True, help="output path prefix (no extension)")
_args = _ap.parse_args()
ROOT = _args.data_dir
OUTDIR = os.path.dirname(os.path.abspath(_args.out_prefix))
os.makedirs(OUTDIR, exist_ok=True)

# ---------- load ----------
region_genes = defaultdict(list)
for row in csv.DictReader(open(os.path.join(ROOT, "all_region_genes.tsv")), delimiter="\t"):
    region_genes[row["region"]].append(dict(
        idx=int(row["gene_idx"]), begin=int(row["begin"]), end=int(row["end"]),
        strand=int(row["strand"]), aa=int(row["aa_len"]),
        gs=min(int(row["begin"]), int(row["end"])), ge=max(int(row["begin"]), int(row["end"]))))

elem_genes = list(csv.DictReader(open(os.path.join(ROOT, "element_gene_meta.tsv")), delimiter="\t"))
for g in elem_genes:
    g["start"] = int(g["start"]); g["end"] = int(g["end"])

links = [r for r in csv.DictReader(open(os.path.join(ROOT, "homology_links.tsv")), delimiter="\t")]

MODULE_COLOR = {
    "lysogeny": "#7B52AB",
    "DNA_packaging": "#2C7FB8",
    "head": "#41B6C4",
    "tail": "#F08A24",
    "lysis": "#D7301F",
    "hypothetical": "#C8C8C8",
    "unassigned": "#E4E4E4",
}

ELEM_REGION = "smbu06_chromosome_1255783-1294647_element_38865bp"

# ---------- module assignment from GenBank product keywords ----------
def keyword_module(product):
    p = product.lower()
    if any(k in p for k in ["terminase", "portal"]):
        return "DNA_packaging"
    if any(k in p for k in ["capsid", "head", "coat", "prohead", "chromosome binding"]):
        return "head"
    if any(k in p for k in ["tail", "baseplate", "tape measure", "tape-measure", "fiber", "structural protein", "base plate", "murein binding"]):
        return "tail"
    if any(k in p for k in ["holin", "amidase", "lysozyme", "lytic", "endolysin", "autolysin", "peptidase"]):
        return "lysis"
    if any(k in p for k in ["integrase", "recombinase", "repressor", "regulator", "sigma", "antitoxin", "toxin", "recombination"]):
        return "lysogeny"
    return "none"


def parse_gb_products(path, region_start):
    """Return list of (rel_start, rel_end, strand, gene, product)."""
    out = []
    cur = None
    prod = None
    gene = None
    for line in open(path):
        m = re.match(r"     CDS             (complement\()?([\d<>]+)\.\.(?>[\d<>]+)\)?", line)
        if m:
            if cur:
                out.append((cur[0], cur[1], cur[2], gene, prod))
            comp = bool(m.group(1))
            endm = re.search(r"\.\.([\d<>]+)", line)
            s = int(m.group(2).replace("<", "").replace(">", "")) - region_start + 1
            e = int(endm.group(1).replace("<", "").replace(">", "")) - region_start + 1
            cur = (min(s, e), max(s, e), -1 if comp else 1)
            prod = None
            gene = None
        elif cur is not None:
            gm = re.search(r'/gene="([^"]+)"', line)
            if gm:
                gene = gm.group(1)
            pm = re.search(r'/product="([^"]+)"', line)
            if pm:
                prod = pm.group(1)
    if cur:
        out.append((cur[0], cur[1], cur[2], gene, prod))
    return out


gb_annot = {
    # NB: efetch with seq_start/seq_stop returns region-relative coordinates (1..N),
    # so region_start is 1 here, not the absolute genome coordinate.
    "Bacillus_subtilis_168_PBSX": parse_gb_products(os.path.join(ROOT, "pbsx_region_Bsub168.gb"), 1),
    "Bacillus_subtilis_168_skin_element_region": parse_gb_products(os.path.join(ROOT, "skin_region_Bsub168.gb"), 1),
}


def module_for_gene(region, g):
    """For PBSX/skin: assign from GenBank product; else None."""
    if region not in gb_annot:
        return None
    best = None
    best_ov = 0
    for (s, e, strand, gene, prod) in gb_annot[region]:
        ov = max(0, min(g["ge"], e) - max(g["gs"], s) + 1)
        if ov > best_ov:
            best_ov = ov
            best = (gene, prod)
    if best and best_ov > 0.5 * (g["ge"] - g["gs"] + 1):
        mod = keyword_module(best[1] or "")
        return mod if mod != "none" else "unassigned", best[0]
    return "unassigned", None


# ---------- orientation + window ----------
def terminase_hit(region):
    for r in links:
        if r["region"] == region and r["query"] == "smbu06GL001272":
            return r
    return None


def flip_of(region):
    h = terminase_hit(region)
    if h is None or region == ELEM_REGION:
        return False
    return int(h["sstart"]) > int(h["send"])


def window_for(region):
    L = max(max(g["begin"], g["end"]) for g in region_genes[region])
    if region == ELEM_REGION:
        return 0, L, False
    if region == "Bacillus_subtilis_168_PBSX":
        return 0, L, False
    if region == "Bacillus_subtilis_168_skin_element_region":
        return 0, L, False
    # homologous span +-2kb in oriented coordinates
    hits = [r for r in links if r["region"] == region]
    flip = flip_of(region)
    PAD = 3000
    if not hits:
        return 0, L, flip
    s0 = min(min(int(r["sstart"]), int(r["send"])) for r in hits)
    s1 = max(max(int(r["sstart"]), int(r["send"])) for r in hits)
    if flip:
        p0, p1 = L - s1 + 1, L - s0 + 1
    else:
        p0, p1 = s0, s1
    return max(0, p0 - PAD), min(L, p1 + PAD), flip


row_specs = [
    (ELEM_REGION, "smbu06 element (this study)", "38,865 bp; 51 CDS (provider annotation)"),
    ("GCF_988228025.1", "E. lactis AF116-27 (GCF_988228025.1)", "terminase 98.6% aa id"),
    ("GCF_049315165.1", "E. lactis BT0171 (GCF_049315165.1)", "terminase 97.3% aa id"),
    ("GCF_049315155.1", "E. durans BT0139 (GCF_049315155.1)", "terminase 93.4% aa id; rev-comp normalised"),
    ("GCF_052655665.1", "E. lactis IS_B39 (GCF_052655665.1)", "terminase 83.7% aa id; rev-comp normalised"),
    ("GCF_006337145.1", "E. faecium HB-1 (GCF_006337145.1)", "terminase 95.0% aa id; rev-comp normalised"),
    ("Enterococcus_faecium_708B_2", "E. faecium 708B_2 (NZ_JARHVH010000002.1)", "draft contig; source of WP_289892796.1"),
    ("Bacillus_subtilis_168_PBSX", "B. subtilis 168 PBSX (NC_000964.3:1,320,570-1,345,795)", "characterised defective prophage (Wood et al. 1990)"),
    ("Bacillus_subtilis_168_skin_element_region", "B. subtilis 168 skin element (NC_000964.3:2,650,900-2,701,000)", "module-architecture comparator; rev-comp normalised"),
]

windows = {region: window_for(region) for region, _, _ in row_specs}

# geometry
FIG_W, FIG_H = 14.0, 8.6
fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
X0 = 0.30  # axes fraction of left edge of sequence area
X1 = 0.94  # right edge
ROW_GAP = 1.0
GENE_H = 0.36

# global scale: bp -> axes fraction, based on the longest displayed window
max_win = max(w[1] - w[0] for w in windows.values())
BP2X = (X1 - X0) / max_win


def tx(region, pos):
    """pos = 1-based coord in original orientation -> axes x"""
    w0, w1, flip = windows[region]
    L = max(max(g["begin"], g["end"]) for g in region_genes[region])
    p = (L - pos + 1) if flip else pos
    return X0 + (p - w0) * BP2X


nrows = len(row_specs)
row_y = {}
for i, (region, _, _) in enumerate(row_specs):
    row_y[region] = nrows - 1 - i

# ---------- draw genes ----------
link_by_region_query = defaultdict(list)
for r in links:
    link_by_region_query[(r["region"], r["query"])].append(r)

# element row colour lookup by pyrodigal gene index
elem_gene_color = {}
for g in region_genes[ELEM_REGION]:
    cx = (g["gs"] + g["ge"]) / 2
    for eg in elem_genes:
        if eg["start"] <= cx <= eg["end"]:
            elem_gene_color[g["idx"]] = MODULE_COLOR.get(eg["module"], MODULE_COLOR["unassigned"])
            break

for region, label, note in row_specs:
    w0, w1, flip = windows[region]
    yy = row_y[region]
    ax.plot([X0, X0 + (w1 - w0) * BP2X], [yy, yy], color="#999999", lw=0.7, zorder=2)
    for g in region_genes[region]:
        if g["ge"] < w0 or g["gs"] > w1:
            continue
        strand = g["strand"] if not flip else -g["strand"]
        x1, x2 = tx(region, g["gs"]), tx(region, g["ge"])
        # apply window clipping
        xa, xb = X0 + (max(w0, min(g["gs"], g["ge"])) - w0) * BP2X if not flip else x1, x2
        xa = max(X0, min(x1, x2))
        xb = min(X0 + (w1 - w0) * BP2X, max(x1, x2))
        w = max(xb - xa, 0.0008)

        if region == ELEM_REGION:
            col = elem_gene_color.get(g["idx"], MODULE_COLOR["unassigned"])
        elif region in gb_annot:
            col, gname = module_for_gene(region, g)
            col = MODULE_COLOR.get(col, MODULE_COLOR["unassigned"]) if col in MODULE_COLOR else MODULE_COLOR["unassigned"]
        else:
            col = MODULE_COLOR["hypothetical"]
            for r in links:
                if r["region"] == region and r["gene_idx"] and int(r["gene_idx"]) == g["idx"] and r["link"] == "1":
                    mod = next((eg["module"] for eg in elem_genes if eg["tag"] == r["query"]), "unassigned")
                    col = MODULE_COLOR.get(mod, MODULE_COLOR["hypothetical"])
                    break

        h = GENE_H
        head = min(w * 0.45, 0.006)
        if strand >= 0:
            pts = [(xa, yy + h / 2), (xb - head, yy + h / 2), (xb, yy), (xb - head, yy - h / 2), (xa, yy - h / 2)]
        else:
            pts = [(xb, yy + h / 2), (xa + head, yy + h / 2), (xa, yy), (xa + head, yy - h / 2), (xb, yy - h / 2)]
        ax.add_patch(Polygon(pts, closed=True, facecolor=col, edgecolor="#404040", lw=0.22, zorder=3))

    ax.text(X0 - 0.008, yy + 0.30, label, ha="right", va="center", fontsize=8.4, fontweight="bold")
    ax.text(X0 - 0.008, yy - 0.30, note, ha="right", va="center", fontsize=6.6, color="#555555", style="italic")

# ---------- links ----------
elem_gene_x = {}
for g in elem_genes:
    cx = (g["start"] + g["end"]) / 2
    elem_gene_x[g["tag"]] = X0 + cx * BP2X

for region, label, note in row_specs:
    if region == ELEM_REGION:
        continue
    yy = row_y[region]
    for r in links:
        if r["region"] != region or r["link"] != "1" or not r["gene_idx"]:
            continue
        gidx = int(r["gene_idx"])
        gg = next((g for g in region_genes[region] if g["idx"] == gidx), None)
        if gg is None:
            continue
        cx = (gg["gs"] + gg["ge"]) / 2
        x2 = tx(region, cx)
        x1 = elem_gene_x.get(r["query"])
        if x1 is None:
            continue
        pid = float(r["pident"])
        lw = 0.30 + (pid - 30) / 70 * 0.75
        alpha = 0.14 + (pid - 30) / 70 * 0.38
        ax.annotate("", xy=(x2, yy + GENE_H / 2), xytext=(x1, row_y[ELEM_REGION] - GENE_H / 2),
                    arrowprops=dict(arrowstyle="-", color="#222222", lw=lw, alpha=alpha,
                                    connectionstyle="arc3,rad=0.0"), zorder=1)

# element gene labels for landmark genes
landmarks = {
    "smbu06GL001245": "int",
    "smbu06GL001248": "rep",
    "smbu06GL001272": "terL",
    "smbu06GL001273": "portal",
    "smbu06GL001284": "tape-measure",
    "smbu06GL001293": "amidase",
}
for g in elem_genes:
    if g["tag"] in landmarks:
        cx = (g["start"] + g["end"]) / 2
        x = X0 + cx * BP2X
        ax.text(x, row_y[ELEM_REGION] + 0.62, landmarks[g["tag"]], ha="center", va="bottom",
                fontsize=6.2, rotation=45, color="#333333")

# population summary on the right for Enterococcus rows
summary = {
    "GCF_988228025.1": "38/51 element proteins matched, 51-100% id",
    "GCF_049315165.1": "38/51 element proteins matched, 51-99% id",
    "GCF_049315155.1": "33/51 element proteins matched, 36-100% id",
    "GCF_052655665.1": "38/51 element proteins matched, 39-100% id",
    "GCF_006337145.1": "39/51 element proteins matched, 51-100% id",
    "Enterococcus_faecium_708B_2": "34/51 element proteins matched, 50-100% id",
    "Bacillus_subtilis_168_PBSX": "no qualifying link; terL 24% id over 97% length",
    "Bacillus_subtilis_168_skin_element_region": "3/51 element proteins matched (yqaJ/yqaK/yqaS, 45-48% id)",
}
for region, txt in summary.items():
    ax.text(X1 + 0.008, row_y[region], txt, ha="left", va="center", fontsize=6.3, color="#333333")

ax.set_xlim(0.0, 1.62)
ax.set_ylim(-1.75, nrows - 0.10)
ax.axis("off")

handles = [mpatches.Patch(facecolor=MODULE_COLOR[k], edgecolor="#404040", label=v)
           for k, v in [("lysogeny", "lysogeny / regulation"), ("DNA_packaging", "DNA packaging"),
                        ("head", "head"), ("tail", "tail / baseplate"), ("lysis", "lysis"),
                        ("hypothetical", "no assignment")]]
leg = ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.62, 0.115), ncol=6,
                fontsize=7.2, frameon=False, title="gene colour = functional module (assigned per row; for PBSX/skin from GenBank annotation keywords)",
                title_fontsize=7.2)
leg._legend_box.align = "left"

sb_len = 5000 * BP2X
ax.plot([X0, X0 + sb_len], [-1.35, -1.35], color="black", lw=1.4)
ax.text(X0 + sb_len / 2, -1.45, "5 kb", ha="center", va="top", fontsize=7.2)

ax.text(X0, -1.05, "Links: real TBLASTN hits (pident ≥ 30%, query coverage ≥ 50%, gene overlap ≥ 50%); grey arrows: no detectable homology to the element.\n"
                   "Colours denote functional module only; they do not imply sequence homology between rows.",
        fontsize=6.4, color="#444444", va="top")

fig.tight_layout()
base = _args.out_prefix
fig.savefig(base + ".pdf")
fig.savefig(base + ".svg")
fig.savefig(base + ".png", dpi=600)
print("saved", base + ".{png,pdf,svg}")
