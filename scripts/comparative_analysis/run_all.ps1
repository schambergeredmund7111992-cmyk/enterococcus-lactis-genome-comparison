# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# run_all.ps1 — 从原始输入重建关键结果（Windows PowerShell）
# 前置：ASCII junction 已建立（mklink /J）：
#   <PROJECT_ROOT>\smbu06 -> 论文3\smbu06
#   <PROJECT_ROOT>\smbu08 -> Downloads\smbu08
#   <PROJECT_ROOT>\comparative_analysis  -> 论文3\smbu06_smbu08_比较基因组完整分析
# 用法: powershell -ExecutionPolicy Bypass -File run_all.ps1 [-SkipDownload] [-SkipMapping] [-SkipTree]

param(
  [switch]$SkipDownload,
  [switch]$SkipMapping,
  [switch]$SkipTree
)

$ErrorActionPreference = 'Continue'
$PY   = '<PYTHON>'
$W    = '<PROJECT_ROOT>\comparative_analysis'
$BASE = '<TOOLS_ROOT>'
$BLAST = "$BASE\ncbi-blast-2.17.0+\bin"
$MM2  = "$BASE\minimap2-win\minimap2-2.31-r1302-windows-x86_64-ucrt64\minimap2.exe"
$MAFFT = "$BASE\mafft-win\mafft.bat"
$IQTREE = "$BASE\iqtree-2.3.6-Windows\bin\iqtree2.exe"

function Step($name, $cmd) {
  Write-Host "=== $name ===" -ForegroundColor Cyan
  $t0 = Get-Date
  Invoke-Expression $cmd
  Write-Host ("--- {0} done in {1:n1} min ---" -f $name, ((Get-Date) - $t0).TotalMinutes) -ForegroundColor Green
}

# 00 manifest
Step '00_inventory'      "$PY `"$W\00_manifest\scripts\inventory.py`""
Step '00_tools_check'    "$PY `"$W\00_manifest\scripts\tools_check.py`""

# 01 assembly QC
Step '01_replicon_stats' "$PY `"$W\01_assembly_qc\scripts\replicon_stats.py`""
if (-not $SkipMapping) {
  Step '01_read_mapping' "$PY `"$W\01_assembly_qc\scripts\read_mapping.py`""
}
Step '01_coverage_analysis' "$PY `"$W\01_assembly_qc\scripts\coverage_analysis.py`""

# 02 taxonomy（refs 需已下载，见 supplementary methods）
Step '02_taxonomy_check' "$PY `"$W\02_taxonomy\scripts\taxonomy_check.py`""

# 03 replicon compare
Step '03_compare_all'    "$PY `"$W\03_replicon_compare\scripts\compare_all.py`""
Step '03_chr_junction'   "$PY `"$W\03_replicon_compare\scripts\chr_junction.py`""
Step '03_pl2_origin'     "$PY `"$W\03_replicon_compare\scripts\pl2_origin.py`""

# 04 plasmidome
Step '04_extract_annotations' "$PY `"$W\04_plasmidome\scripts\extract_annotations.py`""
Step '04_annotate_regions'    "$PY `"$W\04_plasmidome\scripts\annotate_regions.py`""

# 05 bacteriocins
Step '05_build_refs'     "$PY `"$W\05_bacteriocins\scripts\build_refs.py`""
Step '05_nisin_screen'   "$PY `"$W\05_bacteriocins\scripts\nisin_screen.py`""
Step '05_candidate_screen' "$PY `"$W\05_bacteriocins\scripts\candidate_screen.py`""
Step '05_physchem'       "$PY `"$W\05_bacteriocins\scripts\physchem.py`""

# 06 public panel
if (-not $SkipDownload) {
  Step '06_download_panel' "$PY `"$W\06_public_panel\scripts\download_panel_batch.py`""
}
Step '06_panel_analysis' "$PY `"$W\06_public_panel\scripts\panel_analysis.py`""

# 07 phylogenomics
Step '07_core_genes'     "$PY `"$W\07_phylogenomics\scripts\extract_core_genes.py`""
if (-not $SkipTree) {
  Step '07_build_tree'   "$PY `"$W\07_phylogenomics\scripts\build_tree.py`""
}
Step '07_plasmid_network' "$PY `"$W\07_phylogenomics\scripts\plasmid_content_network.py`""

# 08 statistics
Step '08_stats'          "$PY `"$W\08_statistics\scripts\stats_analysis.py`""

# 09 figures（每图独立脚本；QC 汇总）
foreach ($f in @('fig1_overview','fig2_chromosome','fig3_shared_plasmid','fig4_bacteriocin_module',
                 'fig5_plasmid2','figGA_graphical_abstract','figS1_mapping','figS2_nisin','figS3_candidate_alignment')) {
  Step "09_$f" "$PY `"$W\09_figures\scripts\$f.py`""
}
Step '09_fig_qc' "$PY `"$W\09_figures\scripts\fig_qc.py`" $W\09_figures\Fig1_overview\Fig1.pdf $W\09_figures\Fig2_chromosome\Fig2.pdf $W\09_figures\Fig3_shared_plasmid\Fig3.pdf $W\09_figures\Fig4_bacteriocin_module\Fig4.pdf $W\09_figures\Fig5_plasmid2\Fig5.pdf $W\09_figures\GraphicalAbstract\GraphicalAbstract.pdf $W\09_figures\FigS1_mapping\FigS1.pdf $W\09_figures\FigS2_nisin\FigS2.pdf $W\09_figures\FigS3_candidate_alignment\FigS3.pdf"

# 10 supplement
Step '10_supplementary'  "$PY `"$W\10_supplement\scripts\build_supplementary.py`""

Write-Host 'ALL DONE. See 00_manifest/run_log.md and 11_manuscript/ for deliverables.' -ForegroundColor Yellow
