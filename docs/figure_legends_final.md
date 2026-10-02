# Figure legends (v2, PLOS-format style)

**Figure. 1 Overview of the two closed genomes and assembly statistics.** (a, b) Circular representation of the replicons of smbu06 (a; chromosome and Plasmid1) and smbu08 (b; chromosome, Plasmid1 and the additional replicon Plasmid2); highlighted arcs mark the 38,719-bp prophage region (orange; smbu06 chromosome), the excised-element junction attJ (orange; smbu08 chromosome), the bacteriocin-module region (red; shared plasmid) and the two 38,719-bp units of Plasmid2 (purple). (c) Assembly and sequencing statistics recomputed in this study (replicon sizes and GC from FASTA; CDS/tRNA/rRNA from delivered annotations re-parsed; read pairs from FastQC; depths from the provider tables cross-checked by independent mapping, Supplementary Fig. S1).

**Figure. 2 Chromosome comparison: one excisable element and one SNP.** (a) Whole-chromosome dot plot (BLASTN; grey lines = HSPs ≥50 kb, light grey = repeat-family HSPs); the orange bar marks the 38,719-bp segment present only in smbu06 and the arrow the single SNP. (b) Synteny schematic with the prophage block, attL/attR, the smbu08 junction attJ (byte-identical to attR) and the SNP position. (c) The 146-bp junction repeats with the 59 differing positions ticked. (d) Byte-level summary statistics.

**Figure. 3 The 128,837-bp shared plasmid is completely conserved.** (a) Circular map (outer ring: forward-strand CDS; inner ring: reverse-strand CDS; grey curve: GC deviation; red arc: bacteriocin module). (b) Raw-coordinate dot plot of smbu08 versus smbu06 Plasmid1 with the circular-origin break. (c) Rotation-corrected per-kilobase identity (100% across the molecule). (d) Bacteriocin-module neighbourhood. Box: the two independent evidence lines (canonical SHA-256 equality; full-length BLASTN with byte-level verification; 0 SNPs/0 indels) and the absence of a conjugation module.

**Figure. 4 The bacteriocin cassette of the shared plasmid.** (a, b) Gene neighbourhoods of the enterocin P-like locus with the adjacent Sakacin-A-immunity-factor homologue (10-bp intergenic gap) and of the hiracin JM79-like locus. (c) Alignments with characterised homologs (red: differences; boxes: class-IIa motifs). (d) IS-element positions across the shared plasmid; the module is not MGE-enriched relative to the plasmid average. (e) Candidate properties; sequence-level homology only.

**Figure. 5 The smbu08-specific 77.4-kb replicon is the excised prophage.** (a) Model consistent with the observed structures (integrated element in smbu06; excision in smbu08; excised unit persisting as a tandem-dimer replicon). (b) Plasmid2 gene map (102 CDS, modules colour-coded; dashed line = unit boundary with a 1-bp deletion in the second copy). (c) Read depth of the chromosome versus the excised unit (Illumina). (d) Alignment of Plasmid2 to the smbu06 chromosome (reverse strand).

**Figure. 6 Phylogenomic context and distribution in the public panel.** (a) Core-genome maximum-likelihood tree (72 markers; LG+G; SH-aLRT support at major nodes) of the deduplicated panel plus the project genomes, 194 and Lactococcus references, with heatmap columns for plasmid distribution class, cassette-gene presence, nisin status and species. (b) Gene-content relatedness of the shared plasmid (Jaccard distance; average linkage).

**Figure. 7 Functional genome landscape of smbu06.** (a) KEGG pathway classification (top Level-2 categories, coloured by Level-1 class). (b) COG functional categories. (c) GO Level-2 terms (top ten per ontology). smbu08 values are essentially identical (≤2 genes per category); full tables in Supplementary Table S7.

**Figure. 8 Carbohydrate-active enzyme (CAZy) repertoire.** (a) Class counts in both strains (GH = glycoside hydrolases; GT = glycosyltransferases; CBM = carbohydrate-binding modules; CE = carbohydrate esterases; AA = auxiliary activities; PL = polysaccharide lyases). (b) Class distribution (smbu06). (c) Top CAZy families (smbu06). Provider dbCAN-style annotation; no enzyme activity was measured.

**Figure. 9 Safety-related annotation with exact sequence identities.** (a) CARD strict hits (both strains; all four chromosomal): AAC(6′)-Ii 99.45%, msrC 97.15%, eatAv 96.00%, efrA 82.57%. (b) ARDB hits versus the database Min_Identity threshold (all hits below). (c) VFDB category counts (all "Predicted"). Interpretation is limited to sequence-level presence.

**Figure. 10 Pairwise ANIb (%) .** ANIb among the project genomes, strain 194 and comparison genomes, including the strains used in this laboratory's preceding comparative study (E. lactis IDCC 2105 and CX 2-6_2), plus E. faecium 64/3, E. faecalis LD33 and L. lactis MG1363. ANIb: 1,020-bp fragments, ≥70% identity and ≥70% fragment length, bidirectional average; matrix cached in Supplementary Table S7.

**Figure. 11 Pangenome comparison of six genomes.** (a) Gene-level composition of core (in all six genomes), accessory (2–5) and unique (1) orthologue clusters per genome. (b) Pairwise shared orthologue clusters. Proteins predicted uniformly with pyrodigal; clusters built by all-vs-all BLASTP (≥60% identity, ≥60% query coverage) and connected components; comparison panel includes the strains of the preceding study (IDCC 2105, CX 2-6_2).

**Graphical abstract.** Two near-isogenic genome deliveries; a byte-identical 128.8-kb plasmid (two evidence lines) carrying a three-gene bacteriocin cassette; one chromosome-differentiating 38.7-kb PBSX-like element, found excised in smbu08 as a high-copy tandem-dimer replicon; panel-wide distribution shown as counts; cross-host transfer drawn as an explicitly labelled open hypothesis.

## Supplementary figures

**Supplementary Fig. S1 Read assignment, depth and coverage.** (a) Read assignment per dataset. (b, c) Depth profiles of smbu06 (b) and smbu08 (c) replicons. (d) Coverage breadth ≥1× per replicon. (e) Deconvolved depth summary and zero-depth statistics.

**Supplementary Fig. S2 Two-layer nisin search with controls.** Best tblastn identity of the 11 nisin reference proteins against the positive control (F44), negative control (IL1403) and both project genomes; only F44 passes the complete-cluster criterion; the F44 cluster region has zero blastn hits in either genome.

**Supplementary Fig. S3 Candidate peptide alignments.** GL002641 versus enterocin P (94.4%), GL002679 versus hiracin JM79 (59.3%) and the two candidates versus each other; differences in red; class-IIa motifs boxed. Sequence-level only.

**Supplementary Fig. S4 Full annotation of the 38,719-bp element.** All 51 unit CDS with functional-module assignment and NR-derived descriptions (spreadsheet form in Supplementary Table S5).
