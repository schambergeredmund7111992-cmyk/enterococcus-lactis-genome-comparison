# Near-isogenic *Enterococcus lactis* genome comparison (smbu06 vs smbu08)

**Working title**: *Complete genome comparison identifies an excised PBSX-like phage-related element and a conserved bacteriocin-associated plasmid in near-isogenic Enterococcus lactis genomes*

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
* Candidate bacteriocin genes are reported at the **sequence-homology level only**; no expression, processing or activity was measured, and no nisin locus was detected (with positive/negative controls).

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
│   └── longread_validation/      ← independent Nanopore long-read validation
│       ├── 00_manifest/ … 07_report/
├── results/
│   ├── figures/
│   │   ├── main/                 ← Figures 1–11 + GraphicalAbstract (PDF/PNG/SVG)
│   │   ├── supplementary_p3cmp/  ← Supplementary Fig. S1–S4 of the comparative analysis
│   │   ├── supplementary_longread_validation/  ← long-read junction validation figure
│   │   └── FIGURE_INDEX.md       ← maps manuscript figure numbers to files
│   └── tables/
│       ├── supplementary/        ← Supplementary Tables S1–S7 (+ combined workbook)
│       ├── derived_p3cmp/        ← lightweight derived tables per analysis module
│       └── derived_longread_validation/  ← junction counts, depth, model comparison, etc.
└── docs/
    ├── analysis_report_zh.md     ← full Chinese analysis report
    ├── supplementary_methods.md  ← supplementary methods text
    ├── figure_legends_final.md
    ├── claim_audit.tsv / results_to_evidence_map.tsv / reference_verification.md
    └── longread_validation/      ← FINAL_STATUS + validation reports (zh/en)
```

## Software dependencies

Analysis was performed on Windows 11 with Python 3.12.9 (miniconda). Exact versions of all tools are recorded in `software_versions.txt`; the Python packages needed for the scripts in this repository are listed in `requirements.txt` (numpy, pandas, scipy, matplotlib, python-docx, openpyxl).

External command-line tools (not bundled): **minimap2** 2.31, **NCBI BLAST+** 2.17.0, **MAFFT** 7.526, **IQ-TREE** 2.3.6, **NCBI datasets** 18.38, **FastQC**/**MultiQC**, **samtools** 1.24 (optional, for BAM operations; the Windows build used here was obtained from the official MSYS2 `ucrt64` mirror because no official Windows binaries of samtools exist).

## Running the analysis

> **Path placeholders.** The scripts were developed on a local machine and their absolute paths were replaced in this release with placeholders: `<PROJECT_ROOT>` (repository root), `<DATA_ROOT>` (raw-data directory), `<TOOLS_ROOT>` (external bioinformatics tools), `<PYTHON>` (Python interpreter), `<WORKDIR>`. Adapt these at the top of each script before running. Scripts are provided as the exact analysis code; they intentionally reproduce the local directory layout.

1. **Obtain the inputs.** Raw reads and assemblies are not distributed here (see Data availability). Public reference genomes used for taxonomy/phylogeny/panel comparisons are downloaded automatically by `scripts/comparative_analysis/06_public_panel/` (NCBI accessions are recorded in `input_manifest.tsv` and `results/tables/supplementary/TableS3_public_panel/`).
2. **Comparative analysis pipeline.** From `scripts/comparative_analysis/`:
   ```powershell
   powershell -ExecutionPolicy Bypass -File run_all.ps1        # Windows
   sh run_all.sh                                                # Linux / macOS (POSIX companion)
   ```
   `run_all.ps1` / `run_all.sh` document the module order (`00_manifest` → `12_functional_annotation`); switches `-SkipDownload/-SkipMapping/-SkipTree` (PowerShell) or `--skip-download/--skip-mapping/--skip-tree` (POSIX) skip the long-running steps. The analysis scripts themselves are platform-independent Python; only the external binaries and the two runners are platform-specific. Note: on Windows the pipeline expects the ASCII path junctions described in the script header (BLAST+ fails on non-ASCII paths).
3. **Long-read validation.** From `scripts/longread_validation/`, run the scripts in module order:
   `00_manifest/scripts/read_audit2.py` → `01_references/scripts/verify_structure.py` / `build_references.py` → `02_mapping/scripts/map_reads2.py` → `03_supporting_reads/scripts/junction_support2.py` → `per_read_models3.py` → `04_controls/scripts/*` → `05_figures/scripts/fig_s3.py` → `06_tables/scripts/build_tables2.py`.
4. **Expected runtime** (reference machine): comparative pipeline ≈ 2 h end-to-end (dominated by panel download and phylogenetic tree); long-read mapping ≈ 1 h (dominated by SAM output for ~1.2 Gb of Nanopore reads).

## What is *not* included

* Raw sequencing reads (Illumina FASTQ, Nanopore FASTQ/FAST5/POD5) — **not distributed**.
* Complete genome assemblies, GenBank/GFF annotation files of the two project isolates — **not distributed** (unpublished until deposition; see below).
* BAM/BAI/SAM/PAF alignment files, BLAST databases, coverage binary arrays, downloaded public-panel FASTA files and all other files > 25 MB.
* Local execution logs and third-party tool binaries.

## Data availability

> **Raw sequencing data and final assemblies will be deposited in NCBI SRA and GenBank before publication or upon publication.**

Public reference genomes analysed in this repository are available from NCBI under their accessions (see `input_manifest.tsv` and Supplementary Table S3).

## Scope and limitations

All results are computational. Long-read evidence supports the junction sequences and the excised-state minority haplotype but does not distinguish a tandem dimer from a monomer circle; the high copy number of the circular component does not establish autonomous replication; candidate bacteriocin sequences are reported at the homology level; no expression, activity or transfer experiments were performed.

## License

MIT — see `LICENSE`.

## Contact

`TODO: corresponding author name and e-mail` (to be completed before publication)
