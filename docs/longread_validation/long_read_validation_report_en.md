# Long-read validation of the attJ and circular tandem-dimer structure in smbu08

**Date**: 2026-10-02  
**Data**: smbu08 raw Nanopore clean reads (`smbu08.filtered_reads.fq.gz`; 125,722 reads / 1,217.5 Mb / read N50 13,827 bp; SHA-256 `aeafb703…`), smbu08 and smbu06 complete assemblies; smbu06 Nanopore reads used as a positive control.  
**Nature**: computational verification only; no wet-lab work was performed or simulated.

---

## 0. Bottom line

> **Nanopore long reads are consistent with the proposed attJ and circular tandem-dimer structure, but they additionally reveal that the sequenced smbu08 population is a mixture of two chromosome states in which the retained (chr06/194-like, element-present) chromosome predominates (~85%) over the excised attJ state (~15%); the delivered smbu08 chromosome assembly represents the minority (excised) state.**

In short: the attJ state genuinely exists and is strongly supported by model-distinguishing reads (41 strong reads); the plasmid tandem boundary is spanned by ~1,900 long reads. However, the monomer-circle and tandem-dimer interpretations are sequence-equivalent except for a single 1-bp inter-unit indel that the data cannot resolve, and — critically — the sample is a mixed population dominated by the retained chromosome state. This corrects the statement "the element is absent from the smbu08 chromosome": what is element-free is the *assembled haplotype*, not the majority of the population.

## 1. Input audit (`00_manifest/`)

- 125,722 Nanopore reads / 1,217,502,251 bp; mean 9,684 bp; median 6,804 bp; N50 13,827 bp; max 174,754 bp; 48.7% ≥7 kb; 3.5% ≥30 kb. fq.gz SHA-256 `aeafb703ef844c1ba2134e3711161f5fdc734ff35e0c0e9501bc6bccdd7911b9`.
- Assemblies verified by SHA-256 against the original delivery records (chromosome 2,638,187 bp `fb0c3268…`; plasmid FASTA contains 128,837-bp Plasmid1 and 77,437-bp Plasmid2 `9d846394…`).
- Tools: minimap2 2.31-r1302, NCBI BLAST+ 2.17, samtools 1.24 (MSYS2 ucrt64 official mirror build; no official Windows release exists and bioconda has no win-64 build), Python 3.13.
- The raw data had been archived by the user as `<DATA_ROOT>\smbu08.tar.gz`; it was restored from that archive byte-for-byte. **No input file was modified.**

## 2. Structure re-verification (independent second implementation; `01_references/independent_verification.tsv`)

Byte-level, no reliance on aligner boundary placement:

- `chr08 == chr06 with chr06:1,255,783–1,294,501 (38,719 bp) deleted and one base substituted`: left flank 1,255,782 bp with **0 differences**; right flank 1,382,405 bp with exactly **1 SNP** (chr08:2,404,144 C ↔ chr06:2,442,863 T).
- `attL` (chr06 1,255,783–1,255,928), `attJ` (chr08, same coordinates), `attR` (chr06 1,294,502–1,294,647): **attJ is byte-identical to attR**; attL differs from attJ/attR at 59 of 146 positions.
- Element `attL+X+attR` = 38,865 bp (X = 38,573 bp).
- `unit1` (38,719 bp) is **rotation-equivalent to revcomp(attL+X)** (offset 130) — the plasmid unit and the excised segment are exactly the same molecule, not merely similar; `unit2` = unit1 minus 1 bp at unit1 position 36,517 (within a poly-A run).
- All 16 stored reference sequences reproduce byte-for-byte; negative-control provenance verified. The shared 128,837-bp plasmid is rotation-equivalent between assemblies (offset 17,178).

## 3. Junction reference panel

13 references (chr06, chr08, pl1, pl2-as-assembled, attJ window, retained window, unit1/unit2, two 10-kb junction references with ≥3 kb flanks each, two monomer-circle linearisations, a rotated circular-dimer frame, a shuffled negative control and three GC-matched random negative controls), each with construction method, length and SHA-256 (`01_references/reference_construction.tsv`).

## 4. Mapping (`02_mapping/`; commands and logs retained)

minimap2 `-x map-ont` was run twice per reference set (`-c` → PAF; `-a --MD --secondary=yes` → SAM); SAM files were converted to coordinate-sorted, indexed BAM with samtools (1.30 GB / 1.15 GB). Unfiltered outputs, commands, versions, flagstat and idxstats are kept.

- 283,979 alignments: 123,168 primary / 145,097 secondary / 15,714 supplementary; 99.11% of reads mapped.
- **92.2% of primary alignments have mapQ = 0**, because the two chromosomes are identical outside the element and one SNP; chromosome reads split equally between chr06 and chr08. This is a defining property of this near-isogenic dataset and every mapQ-based criterion must state which reference set it refers to.

## 5. Junction-spanning reads (pre-registered criteria)

Criteria fixed before analysis: a single primary alignment continuously spanning the junction (junction covered by match operations), mapQ ≥ 20, matched anchors ≥ 1 kb each side (strict ≥ 3 kb); SAM and PAF parsed independently and cross-checked.

| Junction | Loose reads | Strict reads |
|---|---|---|
| **attJ (chr08, excised state)** | 46 | 26–27 |
| Retained left junction (chr06 flank\|attL) | 257 | 149 |
| Retained right junction (chr06 X\|attR) | 154 | 87 |
| Excised right junction (chr08 attJ\|flank) | 48 | 27 |
| **Plasmid tandem boundary (u1\|u2)** | **1,883** | **992** |
| Plasmid back boundary (origin) | same sequence (both copy boundaries of a circle of identical units share the same junction sequence); covered via the rotated frame: 2,910 spanning alignments |
| GC-matched negative controls (×3) | 248 / 142 / 193 (baseline at ordinary high-depth positions) |
| Shuffled control | 0 |

Strongest attJ reads: `765ab0bd-139d-48a5-a2fa-eb9ba1365e4d` (21,511 bp; anchors 4,984 | 4,985 bp; mapQ 60; identity 99.61%) and `c3ad006f-f202-45ba-b4ec-3bda61199d48` (14,974 bp; anchors 4,987 | 4,979; mapQ 60; identity 99.72%).

## 6. Per-read competing-model analysis (`03_supporting_reads/read_classification3.tsv`)

All 1,678 candidate reads were mapped individually against 10 model references; best alignment scores were compared using the pre-registered thresholds (strong: ΔAS ≥ 10, identity ≥ 90%, coverage ≥ 30%; ambiguous: ΔAS < 10 or identity < 90%).

| Junction class | Strong | Ambiguous | Dominant best model |
|---|---|---|---|
| attJ (chr08) | **41** | 7 | chr08_full / attJ_junction (runner-up chr06_full; ΔAS 3,632–28,336; identity 97–99%) |
| Retained (chr06) | 150 | 108 | chr06_full / retained_window |
| Plasmid forward boundary | **0** | 218 | unit1_to_unit2_junction / pl2 (**ΔAS = 0**) |
| Plasmid back boundary | 0 | 202 | same |
| Negative controls | 0 | 952 | chr08_full (ordinary chromosome positions) |

**Interpretation (to be stated honestly)**: plasmid-boundary reads are abundant, but every one of them ties with the monomer/retained models — because the two tandem copies are identical except for 1 bp, so the plasmid copy boundary, any point of a monomer circle, and the retained chromosome's X region are sequence-equivalent at this level. This is not missing data; it is a *fundamental non-identifiability of the structure from junction-spanning reads*. Discrimination therefore rests on (i) independent existence of the component in the assembly, (ii) the inter-copy 1-bp indel, and (iii) depth.

## 7. Depth, copy number and exclusion of alternatives (`04_controls/`)

Two depth algorithms: samtools depth (primary-only; mapQ thresholds 0 and 20) and an independent CIGAR-walk implementation (primary-only and all-alignments variants). The implementations agree exactly on unique sequence (pl1: 485.3 vs 485.9×); differences in repeat regions are fully explained by secondary-alignment counting semantics.

| Region | primary depth | all alignments | ratio vs chr08 (bootstrap CI95) |
|---|---|---|---|
| chr08 whole | 177.0 | 377.8 | 1.00 |
| chr06 whole | 186.0 | 408.4 | 1.05 (1.04–1.07) |
| **Plasmid2** | **1,002.9** | 2,644.6 | **5.67 (5.51–5.83)** |
| attJ ±2 kb | 73.6 | 153.5 | 0.42 (0.34–0.50) |
| attR ±2 kb (retained-unique) | 241.4 | 697.9 | 1.36 (1.22–1.51) |
| Element X internal | 825.0 | 2,444.0 | 4.66 (4.49–4.83) |
| Plasmid forward boundary ±2 kb | 2,291.9 | 2,792.4 | 12.95 (12.45–13.50) |
| Plasmid origin ±2 kb | 123.4 | 994.4 | 0.70–2.63 (linearisation split, not true low coverage) |
| Negative controls ±2 kb | 148.7–269.9 | 336.2–553.6 | 0.84–1.53 (≈ chromosome level ✓) |

**Copy number (tandem-duplication-corrected)**: weighting each read's alignments by 1/n — element-derived bases 106.5 Mb vs chromosome-derived 907.7 Mb → element abundance ≈ **4.0× chromosome** (CI 3.94–4.05; ≈4 molecule copies per chromosome equivalent under the 77.4-kb dimer interpretation, ≈8 under a 38.7-kb monomer interpretation). **High coverage does not establish autonomous replication.**

Chromosome depth shows a 2.5× ori/ter-like gradient (95×–259×); all window ratios therefore require local normalisation or same-position controls (junction counts are unaffected).

## 8. Read-internal test of the 1-bp inter-unit difference (`04_controls/unit_1bp_*.tsv`)

This indel (unit1 position 36,517, within a poly-A run) is the only sequence-level discriminator between a tandem dimer and a monomer circle at double copy number. Bulk copy1-vs-copy2 distribution comparisons are **non-discriminating** (read assignment between the near-identical copies is decided by alignment score, and the score difference arises from the very indel being tested — a reference-assignment bias; documented in the report). The only clean design is a within-read paired comparison (a single read covering both copies' corresponding windows):

- Exhaustive search across two reference frames (pl2 and the rotated frame): only **2–4 reads** qualify; measured deltas −1, 0, +2, 0.
- Conclusion: **the available data cannot independently validate the 1-bp difference**, and therefore cannot exclude the monomer-circle interpretation on sequence grounds. Its resolution would require assembly-graph information from the provider's pipeline (not available).

## 9. Key finding: smbu08 is a mixed population of two chromosome states (`04_controls/state_junction_counts.tsv`)

Multiple independent lines of evidence:

1. **Same-position junction counts (Nanopore, one read set, identical criteria)**: retained junction (chr06 flank\|attL) **257 loose / 149 strict reads** vs excised junction (chr08 flank\|attJ) **46 loose / 26 strict reads** → retained:excised ≈ **5.6:1 (≈85%:15%)**.
2. **Aligner-independent 20-mer adjacency test** (flank-mer within 300 bp of an attL- or attJ-specific mer on the same read): 184 reads contain the "flank→attL" adjacency (median distance 100 bp); 35 contain "flank→attJ"; the two sets are **mutually exclusive**.
3. **Positive control smbu06 (expected pure-retained)**: 195 vs 8 reads (24:1; ≈4% background) — the method behaves correctly on a pure sample; smbu08's excised signal is far above that background.
4. **Independent technologies and independent prior mappings**: attR (retained-unique) vs attJ (excised-unique) window depths: smbu08 Illumina 196.6 vs 82.3; smbu08 Nanopore 193.3 vs 66.5; smbu06 control attJ only 12.8–29.3 (baseline).
5. **Contamination excluded**: allele counts at the defining SNP (chr08:2,404,144) are C 462 : T 5 (1.07%, at sequencing-error level). If the retained-state molecules came from smbu06/194 contamination (T allele) at ~85% abundance, ~85% T would be observed — it is not. The retained-state molecules carry smbu08's own allele and are internal to this sample.

**Implication**: the delivered smbu08 chromosome assembly (excised state) captured a minority haplotype; the sample majority retains the element. Possible explanations (hypotheses only, pending provider pipeline information): the assembly pipeline resolved the mixed site to the excised path, or a repeat/plasmid-resolution step dropped the prophage from the chromosome.

## 10. Sentences that may / may not be used

**Usable (long-read supported, with qualifiers)**:
> Long-read mappings support the attJ junction state (41 model-distinguishing reads; the retained-state alternative is rejected with ΔAS 3.6–28.3 kb) and place 1,883 reads (≥1 kb anchors, mapQ ≥ 20; 992 strict) across the plasmid's tandem boundary, confirming that the assembled 77.4-kb component is not a junction-level assembly artefact. The tandem-dimer and monomer-circle interpretations are sequence-equivalent except for a single 1-bp inter-unit indel that the available reads cannot resolve.

**Must be added (mixed population)**:
> Long reads show that the sequenced smbu08 population contains both chromosome states, with the retained (element-present) form predominating (257 vs 46 junction-spanning reads, ≈5.6:1); the delivered smbu08 chromosome assembly therefore represents the minority excised haplotype. The retained-state reads carry the smbu08-defining SNP allele, excluding cross-contamination from smbu06/194.

**Not usable (no evidence)**: `induction`, `autonomous replication`, `phage-plasmid`, any equation of high coverage with autonomous replication; `directly validated circular dimer` (monomer vs dimer unresolved); `the element is absent from the smbu08 chromosome` (must become "absent from the delivered assembly / minority haplotype").

## 11. Deliverables and reproducibility

- Figure: `05_figures/{png(600 dpi),pdf,svg}/FigS3_junction_reads.*` (panels A–F; solid = read evidence, dashed = inference).
- Tables: `06_tables/Supplementary_Table_S3_junction_read_support.{xlsx,tsv}`, `junction_reference_manifest.tsv`, `mapping_model_comparison.tsv`, `ambiguous_and_negative_reads.tsv`, `chromosome_state_and_copy_number.tsv`.
- Reads and evidence: `03_supporting_reads/` (candidate and supporting-read FASTA, per-read evidence tables); `02_mapping/bam_subsets/` (±2 kb BAM subsets).
- All scripts under each module's `scripts/`; logs in `00_manifest/logs/`. Every number traces to a command and an output file. Third-party tool provenance in `00_manifest/tool_versions.txt`.

**Limitations**: no official Windows samtools (MSYS2 official-mirror build substituted); no assembly graph available, so the provider's resolution of the mixed site cannot be adjudicated; the 1-bp difference and the monomer/dimer discrimination are statistically underpowered; all conclusions are computational.
