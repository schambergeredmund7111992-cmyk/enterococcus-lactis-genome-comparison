# Figure and table mapping: manuscript (V9) ↔ repository files

Two numbering systems exist in this repository and they **do not coincide**:

* `results/figures/main/Fig1..Fig11` is the historical export of the **complete
  analysis figure set** (11 main figures produced during the analysis, in working
  order). Its numbers are internal and are *not* manuscript numbers.
* The **manuscript** (V9) uses Figure 1–6 and Supplementary Figures S1–S7. Its
  figures were re-laid-out and, for Fig. 2 and Fig. 5, re-drawn; the manuscript
  versions are therefore separate files (see below).

## Manuscript main figures

| Manuscript | Content | Generation script (this repository) | Figure files |
|---|---|---|---|
| Figure 1 | Assembly and sequence architecture of both strains | `scripts/comparative_analysis/09_figures/fig1_overview_revised.py` | working-set export `results/figures/main/Fig1.*` (manuscript layout re-exported in the revision package) |
| Figure 2 | Chromosome-scale structural comparison (38,719-bp differential region; 38,865-bp attL–attR span) | `scripts/manuscript_v9/fig2_chromosome.py` | `results/figures/manuscript_v9/Fig2_manuscript.{png,pdf,svg}` |
| Figure 3 | smbu08 circular component: two units and per-replicon read depth | `scripts/comparative_analysis/09_figures/fig5_plasmid2.py` | working-set export `results/figures/main/Fig5.*` |
| Figure 4 | Shared 128,837-bp plasmid: backbone, loci and identity | `scripts/comparative_analysis/09_figures/fig3_shared_plasmid.py` | working-set export `results/figures/main/Fig2.*` |
| Figure 5 | Candidate bacteriocin-associated loci (structural detail) | `scripts/manuscript_v9/fig5_loci.py` | `results/figures/manuscript_v9/Fig5_manuscript.{png,pdf,svg}` |
| Figure 6 | Public-panel distribution and candidate-gene co-occurrence | `scripts/comparative_analysis/08_statistics/stats_analysis.py` (statistics) and `scripts/comparative_analysis/10_supplement/build_supplementary.py` (table builds); working-set figure export in `results/figures/main/Fig3.*` (lower panel) | see note below |

Note on Figure 6: the manuscript figure combines the plasmid-backbone
distribution (153 of 159 genomes with detectable alignments; categories
near-identical 1 / partial 34 / fragment 88 / trace 30) with the candidate-gene
co-occurrence counts, both computed in `08_statistics`/`06_public_panel`. The
underlying tables are `results/tables/supplementary/TableS3_public_panel/` and
`TableS4_homology_thresholds/`; a stand-alone figure-export script for the exact
manuscript layout is not included in this release.

## Manuscript supplementary figures

| Manuscript | Content | Script | Files |
|---|---|---|---|
| Fig. S1 | Read mapping / contamination summary | `09_figures/figS1_mapping.py` | `results/figures/supplementary_p3cmp/FigS1.*` |
| Fig. S2 | Nisin controls | `09_figures/figS2_nisin.py` | `results/figures/supplementary_p3cmp/FigS2.*` |
| Fig. S3 | Long-read junction support | `scripts/longread_validation/05_figures/` | `results/figures/supplementary_longread_validation/FigS3_junction_reads.*` |
| Fig. S4 | Element annotation (secondary) | `09_figures/figS4_pl2_annotation.py` | `results/figures/supplementary_p3cmp/FigS4.*` |
| Fig. S5 | Pairwise full-length amino-acid alignments of the bacteriocin-like loci | `09_figures/figS3_candidate_alignment.py` | `results/figures/supplementary_p3cmp/FigS3.*` |
| Fig. S6 | Synteny comparison (element vs Enterococcus elements and B. subtilis 168 PBSX/skin) | `scripts/manuscript_v9/plot_synteny.py` | `results/figures/manuscript_v9/FigS6_synteny.*`; inputs in `results/tables/element_tool_validation/synteny/`; caption `FigS6_caption_EN.md` |
| Fig. S7 | ANI heatmap | `09_figures/fig10_ani.py` | `results/figures/manuscript_v9/FigS7_ani.*` (copy of working-set `Fig10.*`) |

## Manuscript supplementary tables

| Manuscript | File |
|---|---|
| S1 | `results/tables/supplementary/TableS1_assembly_stats/` |
| S2 | `results/tables/supplementary/TableS2_candidate_genes/` |
| S3 | `results/tables/supplementary/TableS3_public_panel/` |
| S4 | `results/tables/supplementary/TableS4_homology_thresholds/` |
| S5 | `results/tables/supplementary/TableS5_MGE_and_plasmidome/` |
| S6 | `results/tables/supplementary/TableS6_claim_strength/` |
| S7 | `results/tables/supplementary/TableS7_functional_annotation/` |
| S8 | `results/tables/supplementary/TableS8_ANI_matrix/` |
| S9 | `results/tables/supplementary/TableS9_B168_tblastn/` |
| S10 | `results/tables/supplementary/TableS10_tool_crossvalidation/TableS10_tool_crossvalidation.tsv` |
| S11 | `results/tables/supplementary/TableS11_element_family_members/TableS11_element_family_members.tsv` |
| S12 | `results/tables/supplementary/TableS12_element_span_definitions/TableS12_element_span_definitions.tsv` |
| S13 | `results/tables/supplementary/TableS13_CDS_assignment_audit/TableS13_CDS_assignment_audit.tsv` |

S10–S13 are the manuscript-revision tables; their provenance is described in
`docs/element_tool_validation_methods.md`.
