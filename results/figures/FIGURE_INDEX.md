# Figure index

Mapping between the figure numbers of the current manuscript draft
(`manuscript_fin_修订投稿版`) and the figure files in this repository.
The mapping below was verified by comparing the embedded images in the
manuscript file (by content hash) and their document order with the files
in `results/figures/`.

## Main figures (`results/figures/main/`)

| Manuscript | File stem            | Content |
|------------|----------------------|---------|
| Figure 1   | `Fig1_overview`      | Assembly and read-level overview of smbu06 and smbu08 |
| Figure 2   | `Fig2_chromosome`    | Chromosomal comparison: the 38,719-bp element and one SNP |
| Figure 3   | `Fig5_plasmid2`      | smbu08 circular assembly component (excised tandem-dimer form) |
| Figure 4   | `Fig3_shared_plasmid`| The 128,837-bp plasmid retained in both genomes |
| Figure 5   | `Fig4_bacteriocin_module` | Candidate bacteriocin-associated loci on the shared plasmid |
| Figure 6   | `Fig6_phylogeny`     | Phylogenomic and public-panel context |
| Figure 7   | `Fig7_functional`    | Functional genome landscape (KEGG/COG/GO) |
| Figure 8   | `Fig8_cazy`          | CAZy repertoire |
| Figure 9   | `Fig9_safety`        | Safety-related annotation (CARD/ARDB/VFDB, sequence-level) |
| Figure 10  | `Fig10_ani`          | Pairwise ANIb matrix |
| Figure 11  | `Fig11_pangenome`    | Pangenome comparison |
| —          | `GraphicalAbstract`  | Graphical abstract |

> Note the historical numbering: the figure **file stems** were created earlier
> (before the manuscript was renumbered), so e.g. manuscript Figure 3 is the file
> `Fig5_plasmid2`. The table above is authoritative for the current draft.

Each figure is provided as PDF (vector), PNG (600 dpi) and SVG.

## Supplementary figures

`results/figures/supplementary_p3cmp/` (comparative analysis):

| Provisional | File stem | Content |
|-------------|-----------|---------|
| Fig. S1 | `FigS1` | Read assignment, depth and coverage |
| Fig. S2 | `FigS2` | Two-layer nisin search with positive/negative controls |
| Fig. S3 | `FigS3` | Candidate peptide alignments (enterocin P-like, hiracin JM79-like) |
| Fig. S4 | `FigS4` | Full annotation of the 38,719-bp element (51 CDS) |

`results/figures/supplementary_longread_validation/` (long-read validation):

| Provisional | File stem | Content |
|-------------|-----------|---------|
| Fig. S3 (per current draft) | `FigS3_junction_reads` | Nanopore junction support: structural models and population states (A), junction-spanning reads (B), alignment tracks (C), depth/copy number (D), per-read model competition (E), chromosome-state counts with the smbu06 control (F) |

> **Numbering to confirm before submission:** the current manuscript text refers to
> four supplementary figures (S1 read mapping, S2 nisin controls, S3 long-read
> junction support, S4 secondary annotation). The comparative-analysis set above
> also contains a candidate-alignment figure historically labelled S3. Reconcile
> the final S-numbering before submitting; the file stems are stable.
