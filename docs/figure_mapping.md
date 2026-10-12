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
| Figure 6 | Public-panel distribution and candidate-gene co-occurrence (counts computed from the input tables) | `scripts/manuscript_v9/fig6_panel_distribution.py` | `results/figures/manuscript_v9/Fig6_manuscript.{png,pdf,svg}`; caption `Fig6_caption_EN.md` |

Note on Figure 6: inputs are
`results/tables/supplementary/TableS4_homology_thresholds/pl1_distribution.tsv`
(backbone classes; 153 of 159 panel genomes have detectable alignments:
near-identical 1 / partial 34 / fragment 88 / trace 30) and
`gene_presence_matrix.tsv` (candidate-gene presence at 70%/70% across all 159
genomes; co-occurrence P+I+H 7 / P+I 2 / P only 8 / H only 4). The six panel
genomes without any backbone alignment are the Lactococcus controls F44, G423,
G50, IL1403, IL6288 and IO-1 (verified against
`TableS3_public_panel/panel_metadata.tsv`).

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

> **Historical file names ≠ manuscript numbers.** Several repository directories
> kept their original delivery names (`TableS3_*`, `TableS4_*`, …), which were
> assigned before the manuscript tables were renumbered. Use the mapping below;
> do **not** read a directory name as the manuscript table number.

| Manuscript | Content | Repository file(s) |
|---|---|---|
| S1 | Assembly and read statistics | `results/tables/supplementary/TableS1_assembly_stats/` (replicon stats, coverage per replicon/deconvolved, mapping stats, FastQC/MultiQC summary) |
| S2 | Candidate bacteriocin loci | `results/tables/supplementary/TableS2_candidate_genes/` (alignments of GL002641/GL002679 against the reference precursors, gene table) |
| S3 | **Long-read junction-support evidence** | `results/tables/derived_longread_validation/Supplementary_Table_S3_junction_read_support.tsv` (+ `.xlsx`; derived detail tables alongside: junction counts, model comparison, chromosome-state and copy-number estimates) |
| S4 | **Public-panel metadata, accessions and distributions** | `results/tables/supplementary/TableS3_public_panel/` (accessions, metadata, download records — historical directory name) **plus** `results/tables/supplementary/TableS4_homology_thresholds/` (backbone distribution, gene-presence matrix, co-occurrence/Fisher tables, threshold sensitivity) |
| S5 | Mobile elements and plasmidome annotation | `results/tables/supplementary/TableS5_MGE_and_plasmidome/` |
| S6 | Evidence-strength / homology-threshold audit | `results/tables/supplementary/TableS6_claim_strength/` (chr difference summary, evidence strength, pl2 structure, taxonomy summary) |
| S7 | Functional annotation and safety-related screening | `results/tables/supplementary/TableS7_functional_annotation/` |
| S8 | Complete pairwise ANI matrix | `results/tables/supplementary/TableS8_ANI_matrix/` |
| S9 | B. subtilis 168 TBLASTN comparison of the element proteins | `results/tables/supplementary/TableS9_B168_tblastn/` |
| S10 | Prophage-tool cross-validation | `results/tables/supplementary/TableS10_tool_crossvalidation/TableS10_tool_crossvalidation.tsv` (regenerate with `scripts/element_tool_validation/make_table_s10.py`) |
| S11 | Related elements in public genomes | `results/tables/supplementary/TableS11_element_family_members/TableS11_element_family_members.tsv` |
| S12 | Element span / coordinate definitions | `results/tables/supplementary/TableS12_element_span_definitions/TableS12_element_span_definitions.tsv` |
| S13 | CDS-set assignment audit | `results/tables/supplementary/TableS13_CDS_assignment_audit/TableS13_CDS_assignment_audit.tsv` (regenerate with `scripts/element_tool_validation/make_table_s13_cds_audit.py`) |

Notes:
* `results/tables/supplementary/Supplementary_Tables.xlsx` is the combined
  workbook of the pre-revision tables (S1–S9 content as originally delivered);
  it follows the historical numbering and is **not** authoritative for the
  manuscript table numbers.
* S10–S13 are the manuscript-revision tables; their provenance is described in
  `docs/element_tool_validation_methods.md`.
