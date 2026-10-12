**Supplementary Figure. S6. Gene-order (synteny) comparison of the smbu06 38.7-kb
phage-related element with related enterococcal elements and with the characterised
B. subtilis 168 PBSX and skin elements.**

Rows show, top to bottom: the smbu06 element (1,255,783-1,294,647, the attL..attR span
that contains the 38,719-bp differential region; 51 CDS in the provider annotation);
the terminase-centred regions of five public Enterococcus genomes and of E. faecium
708B_2 (NZ_JARHVH010000002.1), selected from the 159-genome panel by a marker-protein
screen (Supplementary Table S11); the characterised PBSX region of B. subtilis 168
(NC_000964.3:1,320,570-1,345,795; the xkd/xtm region characterised by Wood et al.,
1990) and the B. subtilis 168 skin-element region (NC_000964.3:2,650,900-2,701,000).
Scale bar, 5 kb.

Gene models: open reading frames were predicted de novo with pyrodigal 3.7.1 (meta
mode) for every displayed region, including the element row; the element row's
colouring reflects the provider's 51-CDS module annotation (its label states the
annotation source), not a separate prediction. Arrow direction indicates the coding
strand. Regions whose large terminase is on the minus strand (E. durans BT0139,
E. lactis IS_B39, E. faecium HB-1 and the skin element) are shown reverse-complement
normalised so that every row reads from the lysogeny module towards the lysis module;
this is noted next to each affected row label.

Colours denote functional module only (for the element, from its 51-CDS annotation;
for PBSX/skin, from GenBank /product keywords) and do not imply sequence homology
between rows.

Connecting lines: each line connects one element protein (top row) to its best
homolog in a single reference row and is drawn only where a TBLASTN hit (BLAST+
2.17.0, E <= 1e-5) passes >=30% amino-acid identity, >=50% query coverage and >=50%
overlap with a predicted gene in that row; lines are not pairwise comparisons between
reference rows. The per-row counts at right ("N/51 element proteins matched") count
unique element query proteins with at least one qualifying hit in that row (each
query counted at most once; multiple hits of one query are not summed), with the
identity range of the qualifying hits. Under these thresholds the five Enterococcus
regions and E. faecium 708B_2 match 33-39 of the 51 element proteins at 36-100%
identity, whereas the characterised B. subtilis 168 PBSX region yields no qualifying
link (its large terminase aligns over 97% of the query length at only 24% identity)
and the skin element matches three recombination-associated proteins (YqaJ/YqaK/YqaS,
45-48%). Grey arrows denote genes with no qualifying homology to the element.
Per-hit coordinates, identities and coverages are in the source data
(homology_links.tsv); family accessions and per-protein values are in Supplementary
Table S11.
