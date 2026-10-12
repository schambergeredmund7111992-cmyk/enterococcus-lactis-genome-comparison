# Element classification, tool cross-validation and synteny — methods and status

This document records the revision-stage analyses behind manuscript Supplementary
Tables S10–S13 and Supplementary Figure S6. Everything below was executed as
described; nothing is reported as a result that was not obtained.

## 1. Independent identification (geNomad)

**Software**: geNomad v1.12.0, official marker database v1.9 (2025-02-18 release),
MMseqs2 18.8cc5c. Install via bioconda; database via `genomad download-database`.
**Script**: `scripts/element_tool_validation/run_genomad.py`

Inputs (each run separately; sequences are provider deliveries and are **not**
redistributed here — see `input_manifest.tsv`):
the element sequence (attL..attR, 38,865 bp, chr 1,255,783-1,294,647); the complete
smbu06 and smbu08 assemblies; the 77,437-bp circular component.

```
python run_genomad.py --input <fasta> --out <out_dir> --db <genomad_db> --splits 4 --threads 6
```

Key results (full tables: `results/tables/element_tool_validation/genomad/`;
summary: Table S10):
- smbu06 chromosome — **provirus call** `Chromosome1|provirus_1255929_1301067`
  (45,139 bp; virus score 0.9618; 13 hallmark genes; integrase inside the call).
- smbu08 chromosome — **no call at the corresponding position** (the only
  retained call in either genome is the unrelated 1,036,745-1,047,079 provirus
  present in both). The excision state itself is established by the direct
  sequence comparison of the two delivered assemblies, not by this negative.
- Element alone — virus score 0.9655 (marker-based 0.9965 vs neural-network
  0.5327; both reported); circular component — 0.9703 with 26 hallmarks. All
  classified as Viruses;Duplodnaviria;Heunggongvirae;Uroviricota (phylum);
  Caudoviricetes (class).

**Boundary offsets.** The geNomad 5' provirus boundary coincides with the first
base after attL (1,255,929). The 3' boundary extends **6,420 bp beyond the attR
end of the complete 38,865-bp attL-attR element** (1,294,647) and, equivalently,
**6,566 bp beyond the 3' end of the 38,719-bp differential region** (1,294,501).
Both references are reported; neither is forced to match the att definition.

**Environment note (reproducibility).** The bioconda pairing geNomad 1.7.0 with
MMseqs2 14.7e284 failed reproducibly at the marker-database prefilter step
("Can not allocate entries memory in IndexTable::initMemory"); upgrading to
geNomad 1.12.0 + MMseqs2 18.8cc5c resolved it. The published runs used the
upgraded pair.

## 2. Independent identification (PhiSpy)

**Software**: PhiSpy v5.0.10 (PyPI; Akhter et al. 2012), default settings, without
the optional `--phmms` module.
**Script**: `scripts/element_tool_validation/run_phispy.py`

```
python run_phispy.py --genbank <assembly.gb> --out <out_dir> --threads 4
```

Results (logs and coordinates: `results/tables/element_tool_validation/phispy/`;
summary: Table S10):
- smbu06 — **no retained prophage call**. The highest-ranked candidate
  (Chromosome1 1,264,149-1,279,081; 24 genes, overlapping the element core) was
  dropped with the message "No genes were identified as phage genes". In the
  v5.0.10 evaluation code a candidate is kept only when at least one of its genes
  is classified as a phage gene (default threshold 1); this candidate's count was
  0. The precise contribution of annotation keywords versus sequence features to
  that gene-level classification was **not further determined** and is not
  asserted.
- smbu08 — one retained call, **unrelated to the element** (Plasmid1
  110,154-128,787; IS/recombinase-rich region of the shared plasmid). The
  circular-component regions were high-ranked candidates, dropped for the same
  reason; no candidate at the element position on the chromosome.

The smbu08 negative is reported as tool behaviour only and is not used as
evidence for the excision mechanism.

## 3. Tools not completed (reported as such)

- **PHASTER** — not run in this revision (web-only service). A submission package
  (element ± 10 kb, both assemblies, circular component) accompanies the
  manuscript revision; it contains assembly sequences and is therefore not
  included in this repository. Any returned result (including a negative call)
  belongs in Table S10 as-is.
- **VIBRANT v1.2.1** — installed; the official database setup retrieves Pfam 32.0
  and KEGG 2019-08-10 over `ftp://` endpoints. `ftp://ftp.ebi.ac.uk/...` no longer
  responds and the genome.jp FTP endpoint is non-functional. Reachable HTTPS
  replacements (verified 2026-10-12):
  `https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.hmm.gz`,
  `https://www.genome.jp/ftp/db/kofam/archives/2019-08-10/profiles.tar.gz`;
  the VOG source responds but measured 0-0.17 MB/s from this network. The full
  database (≈ 11 GB final) did not complete. **Not completed — not a negative
  result.**
- **CheckV v1.1.1** — installed; `checkv download_database` measured 0.2-5.5 KB/s
  from the official sources (≈ 1 GB required). **Not completed.** CheckV is a
  quality/completeness assessment against reference genomes, not an independent
  identification vote.

## 4. Related-element survey (Table S11)

**Script**: `scripts/element_tool_validation/panel_family_distribution.py`

Four element marker proteins (integrase GL001245, large terminase GL001272,
portal GL001273, tape-measure GL001284) were searched with TBLASTN (BLAST+
2.17.0, E ≤ 1e-5) against the deduplicated 159-genome panel (accessions in
`results/tables/supplementary/TableS3_public_panel/`; genomes are public and not
redistributed). A genome scored as a related-element carrier when its large
terminase matched at ≥ 60% identity over ≥ 90% of the query.

Result: **19 genomes** (18 public genomes plus the public strain 194 reference
itself), spanning four Enterococcus species (E. lactis, E. durans, E. faecium,
E. raffinosus), large-terminase identity 80.0-98.6% (portal 75.3-96.3%). The
manuscript uses "18 additional public genomes" because the 19th is the
already-public strain 194.

## 5. Synteny comparison (Figure S6)

**Script**: `scripts/manuscript_v9/plot_synteny.py` (inputs:
`results/tables/element_tool_validation/synteny/`).

Reference regions: five Enterococcus element regions (selected by the marker
screen above) + E. faecium 708B_2 (NZ_JARHVH010000002.1:118,000-153,000; the
source genome of WP_289892796.1), the characterised B. subtilis 168 PBSX region
(NC_000964.3:1,320,570-1,345,795; Wood et al. 1990) and the skin-element region
(NC_000964.3:2,650,900-2,701,000). Gene models: pyrodigal 3.7.1 for every
displayed region; the element row's colouring follows the provider 51-CDS
annotation. Homology links require TBLASTN hits (BLAST+ 2.17.0, E ≤ 1e-5) with
≥ 30% identity, ≥ 50% query coverage and ≥ 50% gene overlap; counts are unique
query proteins ("element proteins matched"; each query counted at most once).
Minus-strand regions are reverse-complement normalised for display.
Result: the enterococcal regions share 33-39 of the 51 element proteins at
36-100% identity; the B. subtilis PBSX prototype yields no qualifying link (its
large terminase aligns over 97% of the query at 24% identity).

**Region annotation inputs.** `plot_synteny.py` also reads the two B. subtilis
168 region GenBank files (`pbsx_region_Bsub168.gb`, `skin_region_Bsub168.gb`) to
obtain region sequences and gene models. They are **public** data and are
deliberately not stored in this repository because the repository's general
`*.gb` ignore rule (`.gitignore` line 24) exists to keep unpublished project
assemblies out. Recreate them from NCBI NC_000964.3 with (Entrez direct):

```bash
efetch -db nuccore -id NC_000964.3 -seq_start 1320570 -seq_stop 1345795 \
       -rettype gbwithparts -retmode text > pbsx_region_Bsub168.gb
efetch -db nuccore -id NC_000964.3 -seq_start 2626000 -seq_stop 2701000 \
       -rettype gbwithparts -retmode text > skin_region_Bsub168.gb
```

Verified 2026-10-12 (sha256 prefixes `f31fefde43fa2e26…`, 91,007 bytes;
`71372a957cb4117c…`, 264,184 bytes); the same calls reproduce these files. Place
both files in `results/tables/element_tool_validation/synteny/` (or pass a
`--data-dir` that contains them).

## 6. Element-span and CDS-set definitions (Tables S12, S13)

All coordinates 1-based inclusive, verified against the delivered FASTA
(recomputed 2026-10-12):

| Object | Coordinates | Length |
|---|---|---|
| attL | 1,255,783-1,255,928 | 146 bp |
| internal sequence | 1,255,929-1,294,501 | 38,573 bp |
| attR | 1,294,502-1,294,647 | 146 bp |
| complete attL-attR element | 1,255,783-1,294,647 | 38,865 bp |
| differential region | 1,255,783-1,294,501 | 38,719 bp (= attL + internal) |
| smbu08 attJ | 146 bp, byte-identical to attR | — |
| circular component (Plasmid2) | 77,437 bp = 2 × 38,719 - 1 | — |

The 146-bp difference between the 38,865-bp complete element and the 38,719-bp
differential region is the attR copy retained in smbu08 as attJ.

**CDS-set audit** (`scripts/element_tool_validation/make_table_s13_cds_audit.py`;
output Table S13): 50 of the 51 provider-annotated CDS lie **fully inside** the
internal sequence; the 51st (smbu06GL001295, 114 bp / 38 codons, complement)
lies in the attR region and crosses the 3' end of the complete element.
"the 51 element proteins" denotes this **provider-associated comparison set**
(50 fully internal CDS plus the boundary CDS), not the count of complete CDS
inside the 38,573-bp sequence.

## 7. Classification wording used in the manuscript

"an att-flanked phage-related element with a PBSX-family-annotated terminase".
The terminase similarity (full-length, 441/441 aa, 2 substitutions, E = 0.0, to
the automatically annotated record WP_289892796.1) anchors the terminase-family
annotation; geNomad's assignment is at the viral-class level; the related-element
survey identifies a group, not membership in the classical PBSX group. PBSX-group
membership is stated only as an explicitly limited, family-level inference in the
manuscript Discussion. No induction, infection, particle production, host-DNA
packaging or autonomous replication is claimed.
