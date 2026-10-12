# Figure index

Mapping between the figure numbers of the **current (V9) manuscript** and the
figure files in this repository. The mapping was verified against the rendered
manuscript PDF (figure content read from the embedded images) and the repository
files. See also `docs/figure_mapping.md`.

## Main figures

| Manuscript | Figure files | Generation script |
|------------|--------------|-------------------|
| Figure 1   | working-set export `results/figures/main/Fig1.*` | `scripts/comparative_analysis/09_figures/fig1_overview_revised.py` |
| Figure 2   | `results/figures/manuscript_v9/Fig2_manuscript.*` (redrawn for V9: differential-region vs complete-element labels; "an att-flanked phage-related element") | `scripts/manuscript_v9/fig2_chromosome.py` |
| Figure 3   | working-set export `results/figures/main/Fig5.*` | `scripts/comparative_analysis/09_figures/fig5_plasmid2.py` |
| Figure 4   | working-set export `results/figures/main/Fig2.*` | `scripts/comparative_analysis/09_figures/fig3_shared_plasmid.py` |
| Figure 5   | `results/figures/manuscript_v9/Fig5_manuscript.*` (redesigned for V9: candidate-locus structural detail; replaces the earlier version that duplicated Figure 4's backbone panel) | `scripts/manuscript_v9/fig5_loci.py` |
| Figure 6   | `results/figures/manuscript_v9/Fig6_manuscript.*` (counts computed from the input tables; caption `Fig6_caption_EN.md`) | `scripts/manuscript_v9/fig6_panel_distribution.py` (inputs: `results/tables/supplementary/TableS4_homology_thresholds/{pl1_distribution.tsv,gene_presence_matrix.tsv}`) |

## Supplementary figures

| Manuscript | Figure files | Generation script |
|------------|--------------|-------------------|
| Fig. S1 | `results/figures/supplementary_p3cmp/FigS1.*` | `09_figures/figS1_mapping.py` |
| Fig. S2 | `results/figures/supplementary_p3cmp/FigS2.*` | `09_figures/figS2_nisin.py` |
| Fig. S3 | `results/figures/supplementary_longread_validation/FigS3_junction_reads.*` | `scripts/longread_validation/05_figures/` |
| Fig. S4 | `results/figures/supplementary_p3cmp/FigS4.*` | `09_figures/figS4_pl2_annotation.py` |
| Fig. S5 | `results/figures/supplementary_p3cmp/FigS3.*` (working-set stem) | `09_figures/figS3_candidate_alignment.py` |
| Fig. S6 | `results/figures/manuscript_v9/FigS6_synteny.*` (+ caption `FigS6_caption_EN.md`) | `scripts/manuscript_v9/plot_synteny.py` |
| Fig. S7 | `results/figures/manuscript_v9/FigS7_ani.*` (copy of working-set `Fig10.*`) | `09_figures/fig10_ani.py` |

> **Historical numbering warning.** The working-set exports in
> `results/figures/main/` were created before the manuscript was renumbered, so
> file stems do **not** match manuscript numbers (e.g. manuscript Figure 3 is the
> file `Fig5_plasmid2`). The tables above are authoritative for the current (V9)
> manuscript.

All figures are provided as PDF (vector), PNG and SVG.
