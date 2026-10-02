# Supplementary Methods — exact commands, parameters, versions and database dates

Project: smbu06–smbu08 comparative genome analysis workspace
`smbu06_smbu08_比较基因组完整分析`. All steps are reproducible via `run_all.ps1`.

## 0. Environment

- OS: Windows 11 (10.0.26200). No WSL/Docker/compiler was available; all tools are native Windows builds.
- **ASCII junction staging**: BLAST+ LMDB fails on non-ASCII paths (reproduced:
  `LMDB runtime error`). ASCII junctions were therefore created with
  `mklink /J`: `<PROJECT_ROOT>\{smbu06,smbu08,p3cmp,p3ws}`. All tool I/O uses these ASCII paths.
- Python 3.12.9 (miniconda3); numpy 2.2.5, pandas 2.3.3, scipy 1.15.3, matplotlib 3.10.9,
  python-docx 1.2.0, openpyxl 3.1.5, pyrodigal 3.7.1, MultiQC 1.35.
- NCBI BLAST+ 2.17.0+ (blastn/makeblastdb/blastp/tblastn).
- minimap2 2.31-r1302 (community MSYS2 build, win-ngs/minimap2-windows-build; used because
  bwa/bwa-mem2 have no Windows builds; `-x sr` is the documented equivalent for short-read
  alignment and `-x map-ont` for Nanopore).
- MAFFT 7.526; IQ-TREE 2.3.6; NCBI datasets 18.38.0; FastQC 0.12.1 (Java 21); R 4.6.0.
- Figures: matplotlib, Arial, PNG 600 dpi + PDF + SVG (fonts embedded, `pdf.fonttype=42`).

## 1. Inputs and integrity

- Delivered files enumerated with SHA-256 (`00_manifest/input_inventory.tsv`, 232 files).
  Key hashes (first 32 hex): smbu06.Complete.genome.fasta `adf0b9d4d33aa6bff0ee79de06078886`;
  smbu08.Complete.genome.fasta `21f47d0097a6e4e2e995852ea940aeba`.
- Replicon lengths: smbu06 chromosome 2,676,906 bp / Plasmid1 128,837 bp; smbu08 chromosome
  2,638,187 bp / Plasmid1 128,837 bp / Plasmid2 77,437 bp.
- Company QC (BGI delivery tables): read usage 98.26% (smbu06) / 98.29% (smbu08); reported depths
  chr 460×/440×, Plasmid1 570×/600×, Plasmid2 2,600×; GATK correction 15 indels + 0 SNPs
  (smbu06) and 13 + 0 (smbu08).

## 2. Assembly QC and read mapping (module 01)

- FastQC 0.12.1: `java -Xmx1500m -classpath .;sam-1.103.jar;jbzip2-0.9.jar
  uk.ac.babraham.FastQC.FastQCApplication <fq.gz>`; aggregated with MultiQC 1.35.
  (Note: this FastQC build does not parse `--outdir`/`--threads` on Windows; reports fell back
  to the input directory and were moved into the workspace; input delivery left byte-unchanged.)
- Read statistics: smbu06 4,173,470 pairs, smbu08 4,271,009 pairs (2×150 bp); GC 38%;
  deduplicated fraction 65.7–70.0% (FastQC).
- Mapping: minimap2 `-x sr -c --secondary=no -t 8` (Illumina) and `-x map-ont -c --secondary=no -t 8`
  (Nanopore) against a combined reference of both strains' complete assemblies +
  *L. lactis* 14B4 (GCF_003176835.1) + *L. lactis* MG1363 (GCF_000009425.1).
- Per-base depth computed from PAF CIGAR (M/=/X spans; D not counted), cumulative-sum method.
  Deconvolution: because the two chromosomes/plasmids are near-identical, reads split between the
  two target copies; depths are reported per target and as combined sums
  (`coverage_deconvolved.tsv`).
- Results: 100.0% of reads aligned (all four datasets); Lactococcus-reference hits 0.006%/0.010%
  (Illumina) and 0.05%/0.03% (Nanopore) of reads; combined chromosome depth 443.6×/430.4× (Illumina),
  403.2×/387.6× (Nanopore); Plasmid1 490.1×/517.8×; excised-unit (Plasmid2 + prophage region)
  1,777.6×/1,946.2× (i.e., 4.1–5.0× the chromosome depth).
- Breadth: ≥1× coverage 99.999–100% for all replicons; fraction below ¼ mean depth ≤0.0031%
  (Illumina) and ≤2.49% (Nanopore, Plasmid1).

## 3. Taxonomy check (module 02)

- Reference genomes downloaded with `datasets download genome accession`:
  GCF_056582645.1 (*E. lactis* 194; BioProject PRJNA1433818; Complete Genome; submitter Shenzhen
  University), GCF_003176835.1 (14B4), GCF_000009425.1 (MG1363), GCF_000006865.1 (IL1403).
- Circular-identity tests: SHA-256 on raw and on canonical (minimum rotation over both strands,
  Booth's algorithm). smbu06 chromosome vs 194 chromosome: byte-identical linear representation.
  Both plasmids equal 194's plasmid up to reverse-complement rotation (110,638 bp / 127,816 bp).
- ANIb (JSpeciesWS-style): 1,020-bp fragments, BLASTN, best hit per fragment, identity ≥70% and
  alignment ≥70% of fragment length, bidirectional average. smbu06 vs 194 = 100.00%;
  smbu08 vs 194 = 100.00%; vs *L. lactis* 14B4 = 83.81% / 83.64%.
- 16S: de novo rRNA predictions from the delivery (`*.denovo.rRNA.fasta`); exact comparison and
  BLASTN vs 194; external 16S EF102815 and EF100778 retrieved via NCBI efetch (2026-10-02).
  Both strains' 16S copies are identical to 194 (100%, one copy; 99.81–99.94% across the six
  intragenomic copies); vs EF102815 87.29%, vs EF100778 87.50%.

## 4. Replicon comparison (module 03)

- Canonical circular hashing (above). BLASTN 2.17.0+ megablast for backbone alignment
  (`-evalue 1e-5 -max_hsps 2000`); HSPs ≥50 kb define the colinear backbone, smaller HSPs are
  reported separately as repeat-family hits (`repeat_hsps.tsv`, 1,998 rows).
- Base-level diff: byte comparison over aligned blocks; divergence cores decomposed by
  common prefix/suffix.
- Chromosome results: smbu06 = A(1,255,782) + attL(146) + X′(38,573) + attR(146) + B(1,382,259);
  smbu08 = A + attJ(146, byte-identical to attR) + B. Backbone covers 2,638,187/2,638,187 bp of
  smbu08 (100%); identity 100.000% (back segment) / 99.998% (front); exactly one SNP outside the
  junction: smbu08 2,404,143 C→T (smbu06 2,442,862). attL vs attR differ at 59/146 positions;
  longest common substring 31 bp (5′-AAGAAGCCTTCATGGCCGTTCTGAAAATGGA-3′).
- Plasmid1: raw SHA-256 differ; canonical SHA-256 identical
  (`ce3adc8f8818953d158120d5d6b36a35a2caee3418492049f727b2bb7fc18538`); smbu08 sequence equals
  smbu06 rotated forward by 111,659 bp; BLASTN union covers 128,837/128,837 bp at 100.0%; 0 SNPs.
- Plasmid2: dc-megablast self-alignment: positions 38,720–77,437 match 1–38,719 at 99.997%
  (0 mismatches, 1 gap); single-deletion verification by numpy prefix scan: 1-bp deletion (an A)
  at unit position 36,517 of the second copy. Plasmid2 half1 is a rotation (offset 130) of
  revcomp(smbu06 chromosome 1,255,783–1,294,501); plasmid2 vs smbu06 chromosome:
  38,589 bp at 100.0% (reverse strand) + 38,953 bp at 99.856%.

## 5. Plasmidome / MGE annotation (module 04)

- GBK multi-record parser (custom; CDS + translations; `extract_annotations.py`,
  `annotate_regions.py`). Company NR tables used for best-hit descriptions.
- Prophage region (smbu06 1,255,783–1,294,501): 51 CDS; functional modules assigned by
  keyword rules over product/NR description: lysogeny 4 (site-specific integrase GL001245 99.5%,
  HTH repressor RghR, …), DNA packaging 3 (PBSX-family terminase GL001272 99.5%, portal GL001273
  99.7%), head 3, tail 5 (tape measure GL001284 99.0%, tail/BppU), lysis 3 (holin GL001292 100%,
  XhlA-family hemolysin GL001291 100%, N-acetylmuramoyl-L-alanine amidase GL001293 97.9%),
  hypothetical 33. GC 36.35% vs chromosome 38.47%.
- Plasmid2: 102 CDS = exact 2×51 content (copy 1 = copy 2 by module counts).
- MGE inventory (both strains, all replicons): 192 loci by keyword scan of product + NR
  description (transposase/IS/integrase/recombinase/phage). Plasmid1 is IS-rich (ISS1N/IS6,
  IS3, IS982-ISEfm1, IS431mec, IS257/Tn4003, DDE recombinases); no conjugation module
  (no relaxase/Mob/TraG/MPF detected); pl1-terminal YoeC integrase + IS256 transposase.
- Shared-plasmid module region (31,073–74,688): 14 MGE loci = 0.32/kb vs 40/128.8 kb = 0.31/kb
  plasmid-wide — no IS enrichment around the bacteriocin module.
- Company prophage tooling reported no prophages for either genome (their `Prophage.stat.xls`
  files are empty); the 38.7-kb PBSX-like element was identified here from annotation and
  structure.

## 6. Nisin and bacteriocin mining (module 05)

- Nisin reference set (UniProt, fetched 2026-10-02): nisA V5NV19, nisB Q48673, nisC Q48670,
  nisI Q48671, nisP Q48674, nisR Q07597 (Swiss-Prot; the unreliable Q0GU37 was discarded in the
  earlier project and is not used), nisT Q48669, nisF Q48597, nisE Q48598, nisG Q48599, nisK Q48675.
- Positive control: *L. lactis* F44 nisin cluster region CP024954.1:591,800–606,200 (14,401 bp),
  retrieved 2026-10-02. Negative control: IL1403 (GCF_000006865.1).
- Two-layer search: tblastn (11 proteins, `-evalue 1e-5 -seg no`); blastn (F44 region,
  `-evalue 1e-5`). Cluster call: nisA+nisB+nisC with identity ≥80% and qcov ≥80%.
- Results: F44 positive control = complete cluster (11/11 proteins ≥92.98% identity, 100% qcov);
  IL1403 negative control = not detected; smbu06/smbu08 = **not detected** (only ABC-family
  weak homology: nisF 45.81%, nisR 34.96%, nisT 34.45%, nisK 27.45%; blastn layer zero hits).
- Bacteriocin reference set: UniProt search endpoint, 22 name queries (enterocin, hiracin,
  pediocin, leucocin, sakacin, bavaricin, piscicolin, divercin, carnobacteriocin, lacticin,
  lactococcin, nisin, plantaricin, sublancin, enterolysin, cytolysin, durancin, mundticin,
  garvicin, ubericin, avicin, bacteriocin) + `length:[1 TO 400]`, size=500/query; 3,080 unique
  sequences (fetched 2026-10-02).
- BLASTP: all 2,745 smbu06 CDS translations vs the reference set (`-evalue 1e-3 -seg no`).
- Candidates on Plasmid1 (identical in both strains):
  - GL002641 (36,073–36,237, 54 aa, annotated "Bacteriocin enterocin-P"; UniProt O30434 in GBK
    xref): BLASTP vs enterocin P (O30434) 94.44% identity over full 54 aa (3 differences:
    M1V, YDNGI vs YGNGV). Motif variant YGNGV→YDNGI confirmed at the sequence level.
    Physicochemical (computed): MW 5,767.4 Da; pI 6.74; GRAVY −0.211; net charge −0.17 at pH 7.
  - GL002642 (36,247–36,513, 88 aa; "Sakacin-A immunity factor"), 10-bp intergenic gap
    downstream of GL002641 — immunity-protein candidate (no bacteriocin homology hits).
  - GL002679 (69,485–69,688, 67 aa, "Bacteriocin hiracin-JM79"): vs hiracin JM79 59.26%
    identity over 54 aa (divergent relative); YGNGV motif conserved. MW 7,200.3 Da; pI 9.14;
    GRAVY +0.306; net charge +2.74.
- Short-ORF screen: pyrodigal 3.7.1 (`GeneFinder(meta=False, min_gene=90)`, trained on the
  plasmid sequence) on Plasmid1: 64 ORFs of 30–150 aa; 11 not in the company annotation; only one
  carries a double-glycine motif and none a class-IIa motif → no additional bacteriocin-like
  candidates. Signal-peptide and maturation sites were inferred only by alignment to
  characterised homologs (SignalP unavailable in this environment).

## 7. Public panel (module 06)

- Panel: 213 accessions — 208 from the earlier project panel (accession list archived at
  `06_public_panel/metadata/panel_accessions.txt`) + 5 nisin-positive *L. lactis* controls
  (GCA_000344575.1 IO-1, GCA_000761115.1, GCA_000807375.1, GCA_002804185.1, GCA_002804285.1).
  Download: `datasets download genome accession <25 accs per call> --include genome`
  (2026-10-02); 207/208 retrieved; GCA_056487465.1 failed retries but its RefSeq counterpart
  GCF_056487465.1 (same assembly) was retrieved; dedup removes such pairs.
- Metadata: `datasets summary genome accession <40 accs per call> --as-json-lines`
  → `panel_metadata_all.tsv` (full fields), dedup by GCA/GCF-pair (same assembly) and
  (organism, strain) keeping the highest assembly level → `panel_metadata.tsv`.
- Distribution: BLASTN of the 128,837-bp plasmid against a combined database of all deduped
  panel genomes (contig names prefixed with accession); per-genome query-coverage union and
  length-weighted identity; classes: full-length near-identical (cov ≥95%, id ≥99%),
  full-length divergent (cov ≥95%), partial backbone (≥30%), fragment (≥5%), trace (<5%).
- Gene presence: tblastn of the three candidate proteins vs the same database
  (`-evalue 1e-5 -seg no`); primary thresholds identity ≥70% and query coverage ≥70%;
  sensitivity thresholds 60/70 and 80/80.
- Nisin panel scan: tblastn of the 11 nisin proteins; calls as above.
- Results (deduped panel, N reported in `08_statistics/gene_frequency.tsv`): no genome carries a
  plasmid matching the shared plasmid at ≥80% length; the only full-length match is the 194
  reference itself. Closest relatives share 60.8–67.7% of the backbone at 96.7–98.5% identity;
  34 genomes retain ≥30% of the backbone, 118 only fragments or traces.

## 8. Co-occurrence statistics (module 08)

- Fisher exact tests (SciPy `fisher_exact`, one-sided greater) for candidate-gene pairs at all
  three thresholds; odds ratios with Haldane-Anscombe correction and 95% CI (log method).
  Linkage: same-contig distances from raw tblastn hits. No conclusion is drawn where the number
  of co-carrying genomes is <5.

## 9. Phylogenomics (module 07)

- Marker set: all smbu06 CDS annotated "ribosomal protein" (product or NR description) plus 11
  housekeeping genes (DnaA, DnaK, GroEL, RecA, GyrA, GyrB, RpoA, RpoB, RpoC, EF-Tu, EF-G).
- Extraction: single combined BLAST database (all panel genomes + both project strains + 194 +
  14B4 + MG1363 + IL1403; contig name prefix = accession), one tblastn (`-evalue 1e-10 -seg no`);
  best hit per genome per marker; nucleotide segment extracted, frame-corrected by qstart,
  translated (standard code), internal stops rejected.
- Alignment: MAFFT `--quiet --auto` per marker; columns with >20% gaps removed; markers covering
  ≥80% of genomes retained; genomes present in ≥80% of retained markers retained; concatenated
  supermatrix.
- Tree: IQ-TREE 2.3.6. ML search: `-m LG+G -t <converged tree> -nstop 20 -T 8`. Branch support: SH-aLRT 1000 on the fixed best tree (`-te <tree> -alrt 1000 -T 8`) computed in 48 s; UFBoot was attempted but the full ultrafast-bootstrap run exceeded practical runtime on this 201-taxon x 12,817-site matrix, and SH-aLRT is reported instead (UFBoot is not compatible with a fixed topology).
- Plasmid relatedness: gene-content (Jaccard) distance over shared-plasmid ORF content across
  panel genomes (tblastn presence), average-linkage tree; used instead of a backbone-gene tree
  because no panel genome carries a sufficiently similar full plasmid (max 96.7% identity).

## 10. Figure QC

- All figures exported as PNG 600 dpi + PDF + SVG; programmatic QC with PyMuPDF
  (`fig_qc.py`): text outside page bounds, span-level text overlaps (>55% IoU), font embedding;
  per-figure reports in `09_figures/figure_qc_report.txt`.

## Data and code availability

All scripts, logs, intermediate tables and figures are in
`<PROJECT_ROOT>/smbu06_smbu08_比较基因组完整分析/` and can be re-run with `run_all.ps1`.
Raw sequencing reads and assemblies are the BGI delivery packages (unmodified).
