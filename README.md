# Near-isogenic *Enterococcus lactis* genome comparison (smbu06 vs smbu08)

**Working title**: *Complete genome comparison resolves an att-flanked phage-related element and mixed chromosome states in Enterococcus lactis strains 194 and MGU*

This repository contains the analysis code, figure files and lightweight derived result tables supporting the manuscript above. It is a **reproducible-analysis release**: raw sequencing data and unpublished genome assemblies are not included here (see *Data availability*).

---

## Project aim

Two complete bacterial genome deliveries ("smbu06" and "smbu08", near-isogenic *Enterococcus lactis* isolates) were re-analysed from scratch to resolve their structural differences:

* the 38,719-bp differential chromosomal region and the boundaring `attL`/`attR`/`attJ` attachment sites;
* the 77,437-bp smbu08-specific circular assembly component and its relationship to the excised element;
* the 128,837-bp plasmid shared by both isolates, including a candidate bacteriocin (enterocin P-like / Sakacin-A immunity / hiracin JM79-like) cassette;
* chromosomal single-nucleotide differences, public-panel context and functional annotation;
* an independent **long-read (Nanopore) validation** of the `attJ` junction and the tandem boundary of the circular component.

## Key findings (cautious phrasing)

* **Long-read mappings support the attJ junction state and the unit1→unit2 boundary, but do not distinguish a tandem dimer from a monomer circle.** The two tandem units are identical except for a single 1-bp indel, so junction-spanning reads cannot separate the two interpretations.
* The chromosomes are colinear and differ by one SNP plus the 38,719-bp element; `attJ` is byte-identical to `attR` (59/146 positions differ between `attL` and `attJ`/`attR`).
* The 128,837-bp plasmid is completely conserved between the two deliveries (rotation- and strand-normalised sequence identity, 0 SNPs/indels).
* Long reads additionally indicate that the sequenced smbu08 population contains **both** chromosome states, with the element-retained form predominating (~5.6:1 by junction-spanning read counts); the delivered smbu08 chromosome assembly represents the minority (excised) haplotype. This observation is documented in `docs/longread_validation/`.
* Candidate bacteriocin genes are reported at the **sequence-homology level only**; the two loci are truncated homologues (read-supported in-frame stop codons upstream; no complete precursor ORF), no expression, processing or activity was measured, and no nisin locus was detected (with positive/negative controls).
* The element is described as **"an att-flanked phage-related element with a PBSX-family-annotated terminase"**. The large terminase matches a database-annotated PBSX-family enterococcal terminase at 99.5% amino-acid identity over the full length of both proteins (441/441 aa; WP_289892796.1); geNomad classifies the element region at the viral-class level (*Caudoviricetes*, phylum *Uroviricota*); 18 additional public genomes carry closely related elements (large-terminase identity 80.0–98.6%). PBSX-group membership in a stricter sense is stated only as an explicitly limited, family-level inference. Direct protein-level similarity to the *B. subtilis* 168 PBSX prototype is limited (large terminase 24% identity over 97% of the query; no qualifying synteny link — see `docs/element_tool_validation_methods.md`).

## Repository layout

```
.
├── README.md                     ← this file
├── LICENSE                       ← MIT
├── CITATION.cff                  ← citation metadata (author list to be completed)
├── .gitignore
├── software_versions.txt         ← exact tool versions used
├── requirements.txt              ← Python package dependencies
├── input_manifest.tsv            ← input files / accessions with sizes and SHA-256
├── scripts/
│   ├── comparative_analysis/     ← comparative-genomics pipeline (mirrors the working
│   │   ├── run_all.ps1           ←   directory layout: 00_manifest … 12_functional_annotation)
│   │   └── <module>/             ← module scripts
│   ├── longread_validation/      ← independent Nanopore long-read validation
│   │   ├── 00_manifest/ … 07_report/
│   ├── element_tool_validation/  ← revision-stage independent tool validation
│   │   ├── run_genomad.py        ←   geNomad runs (four inputs)
│   │   ├── run_phispy.py         ←   PhiSpy runs (both genomes)
│   │   ├── make_table_s10.py     ←   aggregates tool outputs into Table S10
│   │   ├── panel_family_distribution.py  ← related-element survey (Table S11)
│   │   └── make_table_s13_cds_audit.py   ← CDS-set audit (Table S13)
│   └── manuscript_v9/            ← manuscript (V9) figure scripts
│       ├── fig2_chromosome.py    ←   Figure 2 (redrawn for V9)
│       ├── fig5_loci.py          ←   Figure 5 (redesigned for V9)
│       └── plot_synteny.py       ←   Figure S6
├── results/
│   ├── figures/
│   │   ├── main/                 ← working-set figure export (Fig1–Fig11 + GraphicalAbstract; internal numbering, NOT manuscript numbers)
│   │   ├── manuscript_v9/        ← manuscript figures added/changed in the revision (Fig2, Fig5, FigS6, FigS7; PDF/PNG/SVG + caption)
│   │   ├── supplementary_p3cmp/  ← Supplementary Fig. S1–S4 of the comparative analysis
│   │   ├── supplementary_longread_validation/  ← long-read junction validation figure
│   │   └── FIGURE_INDEX.md       ← maps manuscript figure numbers to files (see also docs/figure_mapping.md)
│   └── tables/
│       ├── supplementary/        ← Supplementary Tables S1–S13 (+ combined workbook)
│       ├── element_tool_validation/  ← geNomad tables/logs, PhiSpy logs/coordinates, synteny input tables
│       ├── derived_p3cmp/        ← lightweight derived tables per analysis module
│       └── derived_longread_validation/  ← junction counts, depth, model comparison, etc.
└── docs/
    ├── analysis_report_zh.md     ← full Chinese analysis report
    ├── supplementary_methods.md  ← supplementary methods text
    ├── element_tool_validation_methods.md  ← tool validation, synteny and span/CDS definitions (Tables S10–S13, Figure S6)
    ├── figure_mapping.md         ← manuscript figure/table ↔ repository file mapping
    ├── figure_legends_final.md
    ├── claim_audit.tsv / results_to_evidence_map.tsv / reference_verification.md
    └── longread_validation/      ← FINAL_STATUS + validation reports (zh/en)
```

## Software dependencies

Analysis was performed on Windows 11 with Python 3.12.9 (miniconda). Exact versions of all tools are recorded in `software_versions.txt`; the Python packages needed for the scripts in this repository are listed in `requirements.txt` (numpy, pandas, scipy, matplotlib, python-docx, openpyxl).

External command-line tools (not bundled): **minimap2** 2.31, **NCBI BLAST+** 2.17.0, **MAFFT** 7.526, **IQ-TREE** 2.3.6, **NCBI datasets** 18.38, **FastQC**/**MultiQC**, **samtools** 1.24 (optional, for BAM operations; the Windows build used here was obtained from the official MSYS2 `ucrt64` mirror because no official Windows binaries of samtools exist).

## Running the analysis

> **Path placeholders.** The scripts were developed on a local machine and their absolute paths were replaced in this release with placeholders: `<PROJECT_ROOT>` (repository root), `<DATA_ROOT>` (raw-data directory), `<TOOLS_ROOT>` (external bioinformatics tools), `<PYTHON>` (Python interpreter), `<WORKDIR>`. Adapt these at the top of each script before running. Scripts are provided as the exact analysis code; they intentionally reproduce the local directory layout.
>
> **Verification scope.** The revision-stage scripts (`scripts/element_tool_validation/`, `scripts/manuscript_v9/`) were executed in the maintainer's environment against the published inputs, and their outputs reproduce the published tables/figures; small self-tests are provided (`selftest_make_table_s10.py`, `selftest_run_genomad.py`). The main comparative-genomics and long-read pipelines are provided as the exact analysis code and were run in the original working environment, but have not been re-executed as part of this repository update.

1. **Obtain the inputs.** Raw reads and assemblies are not distributed here (see Data availability). Public reference genomes used for taxonomy/phylogeny/panel comparisons are downloaded automatically by `scripts/comparative_analysis/06_public_panel/` (NCBI accessions are recorded in `input_manifest.tsv` and `results/tables/supplementary/TableS3_public_panel/`; historical directory name — this is manuscript Table S4).
2. **Comparative analysis pipeline.** From `scripts/comparative_analysis/`:
   ```powershell
   powershell -ExecutionPolicy Bypass -File run_all.ps1        # Windows
   sh run_all.sh                                                # Linux / macOS (POSIX companion)
   ```
   `run_all.ps1` / `run_all.sh` document the module order (`00_manifest` → `12_functional_annotation`); switches `-SkipDownload/-SkipMapping/-SkipTree` (PowerShell) or `--skip-download/--skip-mapping/--skip-tree` (POSIX) skip the long-running steps. The analysis scripts themselves are platform-independent Python; only the external binaries and the two runners are platform-specific. Note: on Windows the pipeline expects the ASCII path junctions described in the script header (BLAST+ fails on non-ASCII paths).
3. **Long-read validation.** From `scripts/longread_validation/`, run the scripts in module order:
   `00_manifest/scripts/read_audit2.py` → `01_references/scripts/verify_structure.py` / `build_references.py` → `02_mapping/scripts/map_reads2.py` → `03_supporting_reads/scripts/junction_support2.py` → `per_read_models3.py` → `04_controls/scripts/*` → `05_figures/scripts/fig_s3.py` → `06_tables/scripts/build_tables2.py`.
4. **Expected runtime** (reference machine): comparative pipeline ≈ 2 h end-to-end (dominated by panel download and phylogenetic tree); long-read mapping ≈ 1 h (dominated by SAM output for ~1.2 Gb of Nanopore reads).
5. **Independent tool validation (revision stage).** From `scripts/element_tool_validation/` (see `docs/element_tool_validation_methods.md` for versions, commands and the exact status of each tool):
   ```bash
   python run_genomad.py --input <element_or_assembly.fasta> --out <out_dir> --db <genomad_db>   # refuses to delete existing output unless --overwrite is given explicitly
   python run_phispy.py  --genbank <assembly.gb> --out <out_dir>
   python make_table_s10.py --genomad-root <dir> --phispy-root <dir> \
       --target-chromosome Chromosome1 --target-start 1255783 --target-end 1294501 \
       --status-file tool_status_example.tsv --out <TableS10.tsv>
   python panel_family_distribution.py --tblastn <marker_vs_panel.tsv> --metadata <panel_metadata.tsv> --out <TableS11.tsv>
   python make_table_s13_cds_audit.py --element-genes <prophage_genes.tsv> --out <TableS13.tsv>
   python selftest_make_table_s10.py && python selftest_run_genomad.py   # fixture-based checks, no study data needed
   ```
   `make_table_s10.py` parses the actual tool outputs (PhiSpy log candidates/drop messages, prophage coordinates, geNomad summaries) and computes target overlap from the supplied coordinates; study-specific tool status (PHASTER not run; VIBRANT/CheckV not completed) comes from `tool_status_example.tsv`, not from script defaults.
   Tool status for the revision: **geNomad v1.12.0 completed** (provirus call at the element position; no call at the corresponding position of the smbu08 chromosome); **PhiSpy v5.0.10 completed** (element region top-ranked on smbu06 but no retained call — dropped at the gene-level phage-gene filter; one unrelated call on the shared plasmid in smbu08); **PHASTER not run** (submission package prepared for the authors; not included here); **VIBRANT and CheckV not completed** (database sources unavailable from the working network — documented, not negative results).
6. **Manuscript figure scripts.** From `scripts/manuscript_v9/`:
   ```bash
   python fig6_panel_distribution.py \
       --distribution ../../results/tables/supplementary/TableS4_homology_thresholds/pl1_distribution.tsv \
       --presence     ../../results/tables/supplementary/TableS4_homology_thresholds/gene_presence_matrix.tsv \
       --panel-metadata ../../results/tables/supplementary/TableS3_public_panel/panel_metadata.tsv \
       --out-prefix ../../results/figures/manuscript_v9/Fig6_manuscript
   ```
   (All Figure 6 counts are computed from these input tables.) `fig2_chromosome.py` and `plot_synteny.py` take their inputs as arguments; see the script headers and `docs/figure_mapping.md` for the figure ↔ file mapping.
7. **Consistency anchors** (recomputed 2026-10-12 against the delivered FASTA; Tables S12/S13): differential region 38,719 bp (1,255,783–1,294,501); complete attL–attR element 38,865 bp (1,255,783–1,294,647); internal sequence 38,573 bp (1,255,929–1,294,501); attJ/attR 146 bp, byte-identical. The geNomad provirus 3′ boundary extends **6,420 bp** beyond the attR end of the complete element and **6,566 bp** beyond the 3′ end of the differential region. The 51-element-protein set is the provider-associated comparison set (50 CDS fully inside the internal sequence + the boundary CDS GL001295; Table S13).

## What is *not* included

* Raw sequencing reads (Illumina FASTQ, Nanopore FASTQ/FAST5/POD5) — **not distributed**.
* Complete genome assemblies, GenBank/GFF annotation files of the two project isolates — **not distributed** (unpublished until deposition; see below).
* BAM/BAI/SAM/PAF alignment files, BLAST databases, coverage binary arrays, downloaded public-panel FASTA files and all other files > 25 MB.
* Local execution logs of the main comparative/long-read pipelines, and third-party tool binaries. (The revision-stage tool-validation logs **are** included under `results/tables/element_tool_validation/`.)

## Data availability

Raw sequencing reads for the two isolates are deposited in NCBI SRA under the BioProjects listed in the manuscript's Data availability statement (records currently under the NCBI release schedule selected for the submissions); **GenBank accession numbers for the complete assemblies have not yet been assigned**. Public reference genomes analysed in this repository are available from NCBI under their accessions (see `input_manifest.tsv` and Supplementary Table S4). No Zenodo DOI has been assigned yet.

## Scope and limitations

All results are computational. Long-read evidence supports the junction sequences and the excised-state minority haplotype but does not distinguish a tandem dimer from a monomer circle; the high copy number of the circular component does not establish autonomous replication; candidate bacteriocin sequences are reported at the homology level; no expression, activity or transfer experiments were performed.

## License

MIT — see `LICENSE`.

## Contact

`TODO: corresponding author name and e-mail` (to be completed before publication)
