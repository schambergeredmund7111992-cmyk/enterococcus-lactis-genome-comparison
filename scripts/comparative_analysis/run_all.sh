#!/usr/bin/env sh
# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
#
# run_all.sh — POSIX (Linux/macOS) companion to run_all.ps1: reproduces the
# comparative-genomics results step by step. The analysis scripts are pure Python and
# platform-independent; only the external binaries below need platform-specific paths.
#
# Usage:  sh run_all.sh [--skip-download] [--skip-mapping] [--skip-tree]
#
# Prerequisites:
#   - <PROJECT_ROOT>/comparative_analysis : this analysis workspace (mirrored layout)
#   - <DATA_ROOT>/{smbu06,smbu08}         : the delivered data (see input_manifest.tsv)
#   - external tools installed for your platform (BLAST+, minimap2, MAFFT, IQ-TREE, NCBI datasets)

PY='<PYTHON>'
W='<PROJECT_ROOT>/comparative_analysis'
BLAST='<TOOLS_ROOT>/ncbi-blast-2.17.0+/bin'          # makeblastdb, blastn, blastp, tblastn
MM2='<TOOLS_ROOT>/minimap2/bin/minimap2'             # minimap2 (Linux/macOS build)
MAFFT='<TOOLS_ROOT>/mafft/bin/mafft'
IQTREE='<TOOLS_ROOT>/iqtree2/bin/iqtree2'
export PATH="$BLAST:$PATH"

SKIP_DOWNLOAD=0; SKIP_MAPPING=0; SKIP_TREE=0
for a in "$@"; do
  case "$a" in
    --skip-download) SKIP_DOWNLOAD=1 ;;
    --skip-mapping)  SKIP_MAPPING=1 ;;
    --skip-tree)     SKIP_TREE=1 ;;
  esac
done

step() {
  name="$1"; shift
  echo "=== $name ==="
  t0=$(date +%s)
  "$@"
  t1=$(date +%s)
  echo "--- $name done in $(( (t1 - t0) / 60 )) min ---"
}

# 00 manifest
step 00_inventory       "$PY" "$W/00_manifest/scripts/inventory.py"
step 00_tools_check     "$PY" "$W/00_manifest/scripts/tools_check.py"

# 01 assembly QC
step 01_replicon_stats  "$PY" "$W/01_assembly_qc/scripts/replicon_stats.py"
if [ "$SKIP_MAPPING" -eq 0 ]; then
  step 01_read_mapping  "$PY" "$W/01_assembly_qc/scripts/read_mapping.py"
fi
step 01_coverage_analysis "$PY" "$W/01_assembly_qc/scripts/coverage_analysis.py"

# 02 taxonomy (reference downloads required, see supplementary methods)
step 02_taxonomy_check  "$PY" "$W/02_taxonomy/scripts/taxonomy_check.py"

# 03 replicon comparison
step 03_compare_all     "$PY" "$W/03_replicon_compare/scripts/compare_all.py"
step 03_chr_junction    "$PY" "$W/03_replicon_compare/scripts/chr_junction.py"
step 03_pl2_origin      "$PY" "$W/03_replicon_compare/scripts/pl2_origin.py"

# 04 plasmidome
step 04_extract_annotations "$PY" "$W/04_plasmidome/scripts/extract_annotations.py"
step 04_annotate_regions    "$PY" "$W/04_plasmidome/scripts/annotate_regions.py"

# 05 bacteriocins
step 05_build_refs        "$PY" "$W/05_bacteriocins/scripts/build_refs.py"
step 05_nisin_screen      "$PY" "$W/05_bacteriocins/scripts/nisin_screen.py"
step 05_candidate_screen  "$PY" "$W/05_bacteriocins/scripts/candidate_screen.py"
step 05_physchem          "$PY" "$W/05_bacteriocins/scripts/physchem.py"

# 06 public panel
if [ "$SKIP_DOWNLOAD" -eq 0 ]; then
  step 06_download_panel  "$PY" "$W/06_public_panel/scripts/download_panel_batch.py"
fi
step 06_panel_analysis    "$PY" "$W/06_public_panel/scripts/panel_analysis.py"

# 07 phylogenomics
step 07_core_genes        "$PY" "$W/07_phylogenomics/scripts/extract_core_genes.py"
if [ "$SKIP_TREE" -eq 0 ]; then
  step 07_build_tree      "$PY" "$W/07_phylogenomics/scripts/build_tree.py"
fi
step 07_plasmid_network   "$PY" "$W/07_phylogenomics/scripts/plasmid_content_network.py"

# 08 statistics
step 08_stats             "$PY" "$W/08_statistics/scripts/stats_analysis.py"

# 09 figures (one script per figure, then QC summary)
for f in fig1_overview fig2_chromosome fig3_shared_plasmid fig4_bacteriocin_module \
         fig5_plasmid2 figGA_graphical_abstract figS1_mapping figS2_nisin figS3_candidate_alignment; do
  step "09_$f" "$PY" "$W/09_figures/scripts/$f.py"
done
step 09_fig_qc "$PY" "$W/09_figures/scripts/fig_qc.py" \
  "$W/09_figures/Fig1_overview/Fig1.pdf" "$W/09_figures/Fig2_chromosome/Fig2.pdf" \
  "$W/09_figures/Fig3_shared_plasmid/Fig3.pdf" "$W/09_figures/Fig4_bacteriocin_module/Fig4.pdf" \
  "$W/09_figures/Fig5_plasmid2/Fig5.pdf" "$W/09_figures/GraphicalAbstract/GraphicalAbstract.pdf" \
  "$W/09_figures/FigS1_mapping/FigS1.pdf" "$W/09_figures/FigS2_nisin/FigS2.pdf" \
  "$W/09_figures/FigS3_candidate_alignment/FigS3.pdf"

# 10 supplement
step 10_supplementary     "$PY" "$W/10_supplement/scripts/build_supplementary.py"

echo 'ALL DONE. See 00_manifest/run_log.md and 11_manuscript/ for deliverables.'
