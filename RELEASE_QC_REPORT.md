# RELEASE_QC_REPORT — github_release 发布前检查

**检查对象**: `github_release/`（论文投稿用可复现发布目录）  
**检查时间**: 2026-10-02  
**文件总数**: 314 ｜ **总体积**: 28.7 MB

## 1. 通过项

- [通过] 无 > 25 MB 的文件（最大文件 7.7 MB: `results/tables/derived_longread_validation/junction_candidates2.tsv`）
- [通过] 敏感内容扫描（个人绝对路径 / 用户目录名 / 本机 junction 目录名 / AI 助手署名 / token 类 / 邮箱）= 0 处命中
- [通过] README 引用的路径全部存在（核验 8 条，缺失 0 条）
- [通过] 主流程脚本存在且 README 已说明：scripts/comparative_analysis/run_all.ps1（README 提及: True）
- [通过] 已包含：README.md / LICENSE(MIT) / CITATION.cff / .gitignore / requirements.txt / software_versions.txt / input_manifest.tsv
- [通过] 未包含原始测序数据、组装/注释原始文件、BAM/BAI/SAM/PAF、BLAST 数据库、日志、第三方工具链（见第 3 节）

## 2. 待人工确认项

1. **CITATION.cff 与 LICENSE 的作者信息为占位符**（`TODO: add the manuscript author list here`）——投稿前请填写真实作者；本发布不预置任何署名。
2. **README 联系方式为占位符**（`TODO: corresponding author name and e-mail`）。
3. **附图编号**：当前稿件正文按"S1 读段映射 / S2 nisin 对照 / S3 长读段 junction 证据 / S4 次级注释"编号，而比较基因组工作区历史编号中另有 S3（候选肽比对）；投稿前请统一（文件 stem 稳定，见 `results/figures/FIGURE_INDEX.md`）。
4. **脚本路径占位符**：所有脚本中的本机绝对路径已替换为 `<PROJECT_ROOT>`/`<DATA_ROOT>`/`<TOOLS_ROOT>`/`<HOME>` 等占位符（脚本头部有 NOTE 说明），复现前需按本机布局替换。
5. `input_manifest.tsv` 列出的是**未随包分发**的原始输入文件（含 SHA-256），用于核对数据身份；共享原始数据前请确认 SRA/GenBank 已提交。
6. 论文正文/图注中若含未公开的样本来源信息，请在正式公开仓库前再次核对。

## 3. 排除项（未复制/未提交）

- 排除: 原始 Illumina/Nanopore FASTQ/FAST5/POD5（smbu06/smbu08 全部读段数据）
- 排除: 完整组装 FASTA、GenBank（.gb/.gbk）、GFF 等原始注释文件（未公开）
- 排除: BAM/BAI/SAM/PAF 及其子集、BLAST 数据库（*.ndb/nhr/nin/nsq 等）
- 排除: 公共面板基因组下载（panel_all.fna 等）与 07_phylogenomics 中间比对
- 排除: 覆盖率二值数组（*.npz）、MultiQC 完整报告目录中的冗余文件（保留统计表）
- 排除: 本地执行日志（*.log，含个人路径）与 `_superseded_prior_session/` 历史中间文件
- 排除: 第三方工具链（samtools 1.24 MSYS2 包及依赖、bioinfo_tools 等）
- 排除: 论文草稿 docx/manuscript 全文（未纳入本代码发布；如需公开请另行确认）
- 排除: 内部投稿策略笔记（submission_readiness.md 等）与给导师/老师的交付脚本与说明

## 4. 文件清单（全部文件及大小）

| 文件 | 大小 (KB) |
|---|---|
| .gitignore | 1.6 |
| CITATION.cff | 1.5 |
| LICENSE | 1.1 |
| README.md | 8.4 |
| RELEASE_QC_REPORT.md | 27.6 |
| docs/analysis_report_zh.md | 8.6 |
| docs/claim_audit.tsv | 5.0 |
| docs/figure_legends_final.md | 6.6 |
| docs/longread_validation/FINAL_STATUS.md | 7.1 |
| docs/longread_validation/long_read_validation_report.md | 14.6 |
| docs/longread_validation/long_read_validation_report_en.md | 14.1 |
| docs/reference_verification.md | 3.4 |
| docs/results_to_evidence_map.tsv | 3.8 |
| docs/supplementary_methods.md | 14.5 |
| input_manifest.tsv | 40.3 |
| requirements.txt | 0.3 |
| results/figures/FIGURE_INDEX.md | 3.0 |
| results/figures/main/Fig1.pdf | 56.2 |
| results/figures/main/Fig1.png | 518.9 |
| results/figures/main/Fig1.svg | 53.9 |
| results/figures/main/Fig10.pdf | 39.2 |
| results/figures/main/Fig10.png | 535.4 |
| results/figures/main/Fig10.svg | 33.7 |
| results/figures/main/Fig11.pdf | 45.6 |
| results/figures/main/Fig11.png | 458.9 |
| results/figures/main/Fig11.svg | 36.7 |
| results/figures/main/Fig2.pdf | 57.9 |
| results/figures/main/Fig2.png | 527.1 |
| results/figures/main/Fig2.svg | 40.6 |
| results/figures/main/Fig3.pdf | 86.0 |
| results/figures/main/Fig3.png | 829.9 |
| results/figures/main/Fig3.svg | 179.0 |
| results/figures/main/Fig4.pdf | 81.9 |
| results/figures/main/Fig4.png | 626.2 |
| results/figures/main/Fig4.svg | 79.0 |
| results/figures/main/Fig5.pdf | 76.9 |
| results/figures/main/Fig5.png | 535.7 |
| results/figures/main/Fig5.svg | 169.0 |
| results/figures/main/Fig6.pdf | 79.5 |
| results/figures/main/Fig6.png | 1535.6 |
| results/figures/main/Fig6.svg | 473.5 |
| results/figures/main/Fig7.pdf | 43.1 |
| results/figures/main/Fig7.png | 893.8 |
| results/figures/main/Fig7.svg | 63.0 |
| results/figures/main/Fig8.pdf | 44.6 |
| results/figures/main/Fig8.png | 443.8 |
| results/figures/main/Fig8.svg | 37.6 |
| results/figures/main/Fig9.pdf | 48.4 |
| results/figures/main/Fig9.png | 634.1 |
| results/figures/main/Fig9.svg | 33.5 |
| results/figures/main/GraphicalAbstract.pdf | 54.0 |
| results/figures/main/GraphicalAbstract.png | 505.9 |
| results/figures/main/GraphicalAbstract.svg | 14.2 |
| results/figures/supplementary_longread_validation/FigS3_junction_reads.pdf | 60.8 |
| results/figures/supplementary_longread_validation/FigS3_junction_reads.png | 1129.1 |
| results/figures/supplementary_longread_validation/FigS3_junction_reads.svg | 480.0 |
| results/figures/supplementary_p3cmp/FigS1.pdf | 69.7 |
| results/figures/supplementary_p3cmp/FigS1.png | 602.3 |
| results/figures/supplementary_p3cmp/FigS1.svg | 95.4 |
| results/figures/supplementary_p3cmp/FigS2.pdf | 36.4 |
| results/figures/supplementary_p3cmp/FigS2.png | 279.6 |
| results/figures/supplementary_p3cmp/FigS2.svg | 23.2 |
| results/figures/supplementary_p3cmp/FigS3.pdf | 58.9 |
| results/figures/supplementary_p3cmp/FigS3.png | 281.4 |
| results/figures/supplementary_p3cmp/FigS3.svg | 68.9 |
| results/figures/supplementary_p3cmp/FigS4.pdf | 58.8 |
| results/figures/supplementary_p3cmp/FigS4.png | 984.3 |
| results/figures/supplementary_p3cmp/FigS4.svg | 54.2 |
| results/tables/derived_longread_validation/Supplementary_Table_S3_junction_read_support.tsv | 100.4 |
| results/tables/derived_longread_validation/Supplementary_Table_S3_junction_read_support.xlsx | 216.9 |
| results/tables/derived_longread_validation/alignment_evidence_attJ_junction.tsv | 19.0 |
| results/tables/derived_longread_validation/alignment_evidence_neg_ctrl_chr08_1.tsv | 186.3 |
| results/tables/derived_longread_validation/alignment_evidence_neg_ctrl_chr08_2.tsv | 113.0 |
| results/tables/derived_longread_validation/alignment_evidence_neg_ctrl_chr08_3.tsv | 130.5 |
| results/tables/derived_longread_validation/alignment_evidence_retained_window.tsv | 168.5 |
| results/tables/derived_longread_validation/ambiguous_and_negative_reads.tsv | 431.4 |
| results/tables/derived_longread_validation/att_structure_evidence.tsv | 1.6 |
| results/tables/derived_longread_validation/chromosome_state_and_copy_number.tsv | 2.5 |
| results/tables/derived_longread_validation/copy_number_corrected.tsv | 0.6 |
| results/tables/derived_longread_validation/depth_ratio_bootstrap3.tsv | 0.9 |
| results/tables/derived_longread_validation/depth_region_stats3.tsv | 1.3 |
| results/tables/derived_longread_validation/depth_samtools_regions.tsv | 1.4 |
| results/tables/derived_longread_validation/independent_verification.tsv | 4.1 |
| results/tables/derived_longread_validation/input_inventory.tsv | 1.4 |
| results/tables/derived_longread_validation/junction_candidates2.tsv | 7507.0 |
| results/tables/derived_longread_validation/junction_reference_manifest.tsv | 2.6 |
| results/tables/derived_longread_validation/junction_support_summary.tsv | 0.9 |
| results/tables/derived_longread_validation/mapping_model_comparison.tsv | 218.7 |
| results/tables/derived_longread_validation/mixed_structure_counts.tsv | 0.2 |
| results/tables/derived_longread_validation/model_scores2.tsv | 341.8 |
| results/tables/derived_longread_validation/model_scores3.tsv | 620.9 |
| results/tables/derived_longread_validation/multimapping_stats.tsv | 0.3 |
| results/tables/derived_longread_validation/nanopore_read_stats.tsv | 0.4 |
| results/tables/derived_longread_validation/read_classification2.tsv | 135.5 |
| results/tables/derived_longread_validation/read_classification3.tsv | 218.7 |
| results/tables/derived_longread_validation/reference_construction.tsv | 2.6 |
| results/tables/derived_longread_validation/state_junction_counts.tsv | 1.4 |
| results/tables/derived_longread_validation/tool_versions.txt | 0.5 |
| results/tables/derived_longread_validation/unit_1bp_withinread.tsv | 0.5 |
| results/tables/derived_longread_validation/unit_1bp_withinread2.tsv | 0.4 |
| results/tables/derived_longread_validation/unit_1bp_withinread3.tsv | 0.2 |
| results/tables/derived_longread_validation/unit_1bp_withinread_pairs.tsv | 1.3 |
| results/tables/derived_longread_validation/unit_1bp_withinread_pairs2.tsv | 3.2 |
| results/tables/derived_longread_validation/unit_1bp_withinread_pairs3.tsv | 0.3 |
| results/tables/derived_longread_validation/units_and_junctions.tsv | 0.5 |
| results/tables/derived_p3cmp/01_assembly_qc/mapping/contamination_summary.tsv | 0.3 |
| results/tables/derived_p3cmp/01_assembly_qc/mapping/coverage_deconvolved.tsv | 0.3 |
| results/tables/derived_p3cmp/01_assembly_qc/mapping/coverage_per_replicon.tsv | 2.2 |
| results/tables/derived_p3cmp/01_assembly_qc/mapping/coverage_uniformity.tsv | 0.6 |
| results/tables/derived_p3cmp/01_assembly_qc/mapping/mapping_stats.tsv | 0.5 |
| results/tables/derived_p3cmp/01_assembly_qc/reads_qc/multiqc_report/multiqc_data/fastqc-status-check-heatmap.txt | 0.4 |
| results/tables/derived_p3cmp/01_assembly_qc/reads_qc/multiqc_report/multiqc_data/fastqc_per_base_n_content_plot.txt | 3.3 |
| results/tables/derived_p3cmp/01_assembly_qc/reads_qc/multiqc_report/multiqc_data/fastqc_per_base_sequence_quality_plot.txt | 3.8 |
| results/tables/derived_p3cmp/01_assembly_qc/reads_qc/multiqc_report/multiqc_data/fastqc_per_sequence_gc_content_plot_Counts.txt | 5.9 |
| results/tables/derived_p3cmp/01_assembly_qc/reads_qc/multiqc_report/multiqc_data/fastqc_per_sequence_gc_content_plot_Percentages.txt | 10.2 |
| results/tables/derived_p3cmp/01_assembly_qc/reads_qc/multiqc_report/multiqc_data/fastqc_per_sequence_quality_scores_plot.txt | 1.2 |
| results/tables/derived_p3cmp/01_assembly_qc/reads_qc/multiqc_report/multiqc_data/fastqc_sequence_counts_plot.txt | 0.2 |
| results/tables/derived_p3cmp/01_assembly_qc/reads_qc/multiqc_report/multiqc_data/fastqc_sequence_duplication_levels_plot.txt | 1.5 |
| results/tables/derived_p3cmp/01_assembly_qc/reads_qc/multiqc_report/multiqc_data/llms-full.txt | 40.1 |
| results/tables/derived_p3cmp/01_assembly_qc/reads_qc/multiqc_report/multiqc_data/multiqc_citations.txt | 0.1 |
| results/tables/derived_p3cmp/01_assembly_qc/reads_qc/multiqc_report/multiqc_data/multiqc_fastqc.txt | 1.2 |
| results/tables/derived_p3cmp/01_assembly_qc/reads_qc/multiqc_report/multiqc_data/multiqc_general_stats.txt | 0.4 |
| results/tables/derived_p3cmp/01_assembly_qc/reads_qc/multiqc_report/multiqc_data/multiqc_software_versions.txt | 0.0 |
| results/tables/derived_p3cmp/01_assembly_qc/reads_qc/multiqc_report/multiqc_data/multiqc_sources.txt | 0.6 |
| results/tables/derived_p3cmp/01_assembly_qc/replicon_stats.tsv | 0.6 |
| results/tables/derived_p3cmp/02_taxonomy/16S_comparison.tsv | 0.7 |
| results/tables/derived_p3cmp/02_taxonomy/16S_vs_194.tsv | 1.5 |
| results/tables/derived_p3cmp/02_taxonomy/ani_results.tsv | 0.2 |
| results/tables/derived_p3cmp/02_taxonomy/taxonomy_summary.tsv | 0.6 |
| results/tables/derived_p3cmp/03_replicon_compare/att_region_analysis.tsv | 1.0 |
| results/tables/derived_p3cmp/03_replicon_compare/chr_difference_summary.tsv | 0.7 |
| results/tables/derived_p3cmp/03_replicon_compare/chr_snp_list.tsv | 0.1 |
| results/tables/derived_p3cmp/03_replicon_compare/chromosome_comparison_summary.tsv | 0.3 |
| results/tables/derived_p3cmp/03_replicon_compare/difference_regions.tsv | 0.2 |
| results/tables/derived_p3cmp/03_replicon_compare/evidence_strength.tsv | 1.7 |
| results/tables/derived_p3cmp/03_replicon_compare/identity_hashes.tsv | 0.8 |
| results/tables/derived_p3cmp/03_replicon_compare/pl2_structure_summary.tsv | 0.5 |
| results/tables/derived_p3cmp/03_replicon_compare/pl2_vs_chr06.hsp.tsv | 2.1 |
| results/tables/derived_p3cmp/03_replicon_compare/plasmid_comparison_summary.tsv | 0.5 |
| results/tables/derived_p3cmp/03_replicon_compare/repeat_hsps.tsv | 118.3 |
| results/tables/derived_p3cmp/03_replicon_compare/snp_list.tsv | 0.1 |
| results/tables/derived_p3cmp/04_plasmidome/annotation/full_gene_table.tsv | 546.5 |
| results/tables/derived_p3cmp/04_plasmidome/annotation/mge_inventory.tsv | 26.5 |
| results/tables/derived_p3cmp/04_plasmidome/annotation/phage_module_table.tsv | 9.7 |
| results/tables/derived_p3cmp/04_plasmidome/annotation/pl1_genes.tsv | 26.8 |
| results/tables/derived_p3cmp/04_plasmidome/annotation/pl1_key_genes.tsv | 23.9 |
| results/tables/derived_p3cmp/04_plasmidome/annotation/pl2_annotation_full.tsv | 17.0 |
| results/tables/derived_p3cmp/04_plasmidome/annotation/pl2_genes.tsv | 8.2 |
| results/tables/derived_p3cmp/04_plasmidome/annotation/prophage_annotation_full.tsv | 9.7 |
| results/tables/derived_p3cmp/04_plasmidome/annotation/prophage_genes.tsv | 5.5 |
| results/tables/derived_p3cmp/04_plasmidome/annotation/replicon_gene_summary.tsv | 0.2 |
| results/tables/derived_p3cmp/05_bacteriocins/refs_manifest.tsv | 0.6 |
| results/tables/derived_p3cmp/05_bacteriocins/search/alignment_GL002641_vs_enterocinP.txt | 0.2 |
| results/tables/derived_p3cmp/05_bacteriocins/search/alignment_GL002679_vs_hiracinJM79.txt | 0.2 |
| results/tables/derived_p3cmp/05_bacteriocins/search/candidate_table.tsv | 0.4 |
| results/tables/derived_p3cmp/05_bacteriocins/search/mge_distance_table.tsv | 8.1 |
| results/tables/derived_p3cmp/05_bacteriocins/search/neighborhood_GL002641.tsv | 2.5 |
| results/tables/derived_p3cmp/05_bacteriocins/search/neighborhood_GL002679.tsv | 2.5 |
| results/tables/derived_p3cmp/05_bacteriocins/search/nisin_control_validation.tsv | 0.2 |
| results/tables/derived_p3cmp/05_bacteriocins/search/nisin_status_matrix.tsv | 0.6 |
| results/tables/derived_p3cmp/05_bacteriocins/search/nisin_tblastn_hits.tsv | 32.8 |
| results/tables/derived_p3cmp/05_bacteriocins/search/physchem_candidates.tsv | 0.7 |
| results/tables/derived_p3cmp/05_bacteriocins/search/pl1_vs_bacteriocin_refset.tsv | 132.7 |
| results/tables/derived_p3cmp/05_bacteriocins/search/small_orf_screen.tsv | 6.0 |
| results/tables/derived_p3cmp/06_public_panel/analysis/gene_hits_all.tsv | 7.8 |
| results/tables/derived_p3cmp/06_public_panel/analysis/gene_presence_matrix.tsv | 21.7 |
| results/tables/derived_p3cmp/06_public_panel/analysis/nisin_panel_scan.tsv | 29.4 |
| results/tables/derived_p3cmp/06_public_panel/analysis/nisin_panel_status.tsv | 5.0 |
| results/tables/derived_p3cmp/06_public_panel/analysis/pl1_distribution.tsv | 12.0 |
| results/tables/derived_p3cmp/06_public_panel/metadata/download_log.tsv | 30.2 |
| results/tables/derived_p3cmp/06_public_panel/metadata/download_record.tsv | 18.3 |
| results/tables/derived_p3cmp/06_public_panel/metadata/panel_accessions.txt | 3.5 |
| results/tables/derived_p3cmp/06_public_panel/metadata/panel_metadata.tsv | 31.5 |
| results/tables/derived_p3cmp/06_public_panel/metadata/panel_metadata_all.tsv | 41.0 |
| results/tables/derived_p3cmp/07_phylogenomics/core_gene_table.tsv | 907.5 |
| results/tables/derived_p3cmp/07_phylogenomics/trees/core_tree.treefile | 6.7 |
| results/tables/derived_p3cmp/07_phylogenomics/trees/core_tree_alrt.treefile | 9.0 |
| results/tables/derived_p3cmp/07_phylogenomics/trees/plasmid_gene_content_matrix.tsv | 4.3 |
| results/tables/derived_p3cmp/08_statistics/cooccurrence_fisher.tsv | 1.0 |
| results/tables/derived_p3cmp/08_statistics/gene_frequency.tsv | 1.0 |
| results/tables/derived_p3cmp/08_statistics/linkage_641_679.tsv | 2.1 |
| results/tables/derived_p3cmp/08_statistics/plasmid_gene_association.tsv | 0.7 |
| results/tables/derived_p3cmp/12_functional_annotation/summary/ani_matrix.tsv | 1.3 |
| results/tables/derived_p3cmp/12_functional_annotation/summary/ardb_hits.tsv | 2.9 |
| results/tables/derived_p3cmp/12_functional_annotation/summary/card_hits.tsv | 2.0 |
| results/tables/derived_p3cmp/12_functional_annotation/summary/cazy_5class.tsv | 0.1 |
| results/tables/derived_p3cmp/12_functional_annotation/summary/cazy_families.tsv | 1.8 |
| results/tables/derived_p3cmp/12_functional_annotation/summary/cog.tsv | 2.7 |
| results/tables/derived_p3cmp/12_functional_annotation/summary/go_level2.tsv | 3.0 |
| results/tables/derived_p3cmp/12_functional_annotation/summary/kegg_level1.tsv | 0.4 |
| results/tables/derived_p3cmp/12_functional_annotation/summary/kegg_level2.tsv | 4.3 |
| results/tables/derived_p3cmp/12_functional_annotation/summary/pangenome_summary.tsv | 1.0 |
| results/tables/derived_p3cmp/12_functional_annotation/summary/t3ss_summary.tsv | 0.1 |
| results/tables/derived_p3cmp/12_functional_annotation/summary/vfdb_category_counts.tsv | 0.6 |
| results/tables/derived_p3cmp/12_functional_annotation/summary/vfdb_hits.tsv | 22.6 |
| results/tables/supplementary/Supplementary_Tables.xlsx | 177.0 |
| results/tables/supplementary/TableS1_assembly_stats/contamination_summary.tsv | 0.3 |
| results/tables/supplementary/TableS1_assembly_stats/coverage_deconvolved.tsv | 0.3 |
| results/tables/supplementary/TableS1_assembly_stats/coverage_per_replicon.tsv | 2.2 |
| results/tables/supplementary/TableS1_assembly_stats/coverage_uniformity.tsv | 0.6 |
| results/tables/supplementary/TableS1_assembly_stats/mapping_stats.tsv | 0.5 |
| results/tables/supplementary/TableS1_assembly_stats/multiqc_fastqc.txt | 1.2 |
| results/tables/supplementary/TableS1_assembly_stats/replicon_stats.tsv | 0.6 |
| results/tables/supplementary/TableS2_candidate_genes/alignment_GL002641_vs_enterocinP.txt | 0.2 |
| results/tables/supplementary/TableS2_candidate_genes/alignment_GL002679_vs_hiracinJM79.txt | 0.2 |
| results/tables/supplementary/TableS2_candidate_genes/candidate_table.tsv | 0.4 |
| results/tables/supplementary/TableS2_candidate_genes/neighborhood_GL002641.tsv | 2.5 |
| results/tables/supplementary/TableS2_candidate_genes/neighborhood_GL002679.tsv | 2.5 |
| results/tables/supplementary/TableS2_candidate_genes/physchem_candidates.tsv | 0.7 |
| results/tables/supplementary/TableS3_public_panel/download_log.tsv | 30.2 |
| results/tables/supplementary/TableS3_public_panel/download_record.tsv | 18.3 |
| results/tables/supplementary/TableS3_public_panel/panel_metadata.tsv | 31.5 |
| results/tables/supplementary/TableS3_public_panel/panel_metadata_all.tsv | 41.0 |
| results/tables/supplementary/TableS4_homology_thresholds/cooccurrence_fisher.tsv | 1.0 |
| results/tables/supplementary/TableS4_homology_thresholds/gene_frequency.tsv | 1.0 |
| results/tables/supplementary/TableS4_homology_thresholds/gene_presence_matrix.tsv | 21.7 |
| results/tables/supplementary/TableS4_homology_thresholds/linkage_641_679.tsv | 2.1 |
| results/tables/supplementary/TableS4_homology_thresholds/nisin_panel_scan.tsv | 29.4 |
| results/tables/supplementary/TableS4_homology_thresholds/nisin_panel_status.tsv | 5.0 |
| results/tables/supplementary/TableS4_homology_thresholds/pl1_distribution.tsv | 12.0 |
| results/tables/supplementary/TableS4_homology_thresholds/plasmid_gene_association.tsv | 0.7 |
| results/tables/supplementary/TableS5_MGE_and_plasmidome/mge_distance_table.tsv | 8.1 |
| results/tables/supplementary/TableS5_MGE_and_plasmidome/mge_inventory.tsv | 26.5 |
| results/tables/supplementary/TableS5_MGE_and_plasmidome/pl1_key_genes.tsv | 23.9 |
| results/tables/supplementary/TableS5_MGE_and_plasmidome/pl2_annotation_full.tsv | 17.0 |
| results/tables/supplementary/TableS5_MGE_and_plasmidome/pl2_genes.tsv | 8.2 |
| results/tables/supplementary/TableS5_MGE_and_plasmidome/prophage_annotation_full.tsv | 9.7 |
| results/tables/supplementary/TableS5_MGE_and_plasmidome/prophage_genes.tsv | 5.5 |
| results/tables/supplementary/TableS6_claim_strength/chr_difference_summary.tsv | 0.7 |
| results/tables/supplementary/TableS6_claim_strength/evidence_strength.tsv | 1.7 |
| results/tables/supplementary/TableS6_claim_strength/pl2_structure_summary.tsv | 0.5 |
| results/tables/supplementary/TableS6_claim_strength/taxonomy_summary.tsv | 0.6 |
| results/tables/supplementary/TableS7_functional_annotation/ani_matrix.tsv | 1.3 |
| results/tables/supplementary/TableS7_functional_annotation/ardb_hits.tsv | 2.9 |
| results/tables/supplementary/TableS7_functional_annotation/card_hits.tsv | 2.0 |
| results/tables/supplementary/TableS7_functional_annotation/cazy_5class.tsv | 0.1 |
| results/tables/supplementary/TableS7_functional_annotation/cazy_families.tsv | 1.8 |
| results/tables/supplementary/TableS7_functional_annotation/cog.tsv | 2.7 |
| results/tables/supplementary/TableS7_functional_annotation/go_level2.tsv | 3.0 |
| results/tables/supplementary/TableS7_functional_annotation/kegg_level1.tsv | 0.4 |
| results/tables/supplementary/TableS7_functional_annotation/kegg_level2.tsv | 4.3 |
| results/tables/supplementary/TableS7_functional_annotation/pangenome_summary.tsv | 1.0 |
| results/tables/supplementary/TableS7_functional_annotation/vfdb_category_counts.tsv | 0.6 |
| scripts/comparative_analysis/00_manifest/inventory.py | 4.6 |
| scripts/comparative_analysis/00_manifest/p3lib.py | 6.4 |
| scripts/comparative_analysis/00_manifest/tools_check.py | 3.1 |
| scripts/comparative_analysis/01_assembly_qc/coverage_analysis.py | 4.7 |
| scripts/comparative_analysis/01_assembly_qc/read_mapping.py | 9.6 |
| scripts/comparative_analysis/01_assembly_qc/replicon_stats.py | 6.7 |
| scripts/comparative_analysis/02_taxonomy/taxonomy_check.py | 14.7 |
| scripts/comparative_analysis/03_replicon_compare/chr_junction.py | 5.2 |
| scripts/comparative_analysis/03_replicon_compare/compare_all.py | 13.9 |
| scripts/comparative_analysis/03_replicon_compare/pl2_origin.py | 3.4 |
| scripts/comparative_analysis/04_plasmidome/annotate_regions.py | 8.1 |
| scripts/comparative_analysis/04_plasmidome/extract_annotations.py | 7.8 |
| scripts/comparative_analysis/05_bacteriocins/build_refs.py | 5.5 |
| scripts/comparative_analysis/05_bacteriocins/candidate_screen.py | 12.5 |
| scripts/comparative_analysis/05_bacteriocins/nisin_screen.py | 6.0 |
| scripts/comparative_analysis/05_bacteriocins/physchem.py | 5.7 |
| scripts/comparative_analysis/06_public_panel/download_panel.py | 3.5 |
| scripts/comparative_analysis/06_public_panel/download_panel_batch.py | 3.9 |
| scripts/comparative_analysis/06_public_panel/panel_analysis.py | 13.5 |
| scripts/comparative_analysis/06_public_panel/retry_missing.py | 2.5 |
| scripts/comparative_analysis/07_phylogenomics/build_tree.py | 4.3 |
| scripts/comparative_analysis/07_phylogenomics/extract_core_genes.py | 6.9 |
| scripts/comparative_analysis/07_phylogenomics/plasmid_content_network.py | 4.4 |
| scripts/comparative_analysis/08_statistics/stats_analysis.py | 6.6 |
| scripts/comparative_analysis/09_figures/fig10_ani.py | 5.8 |
| scripts/comparative_analysis/09_figures/fig11_pangenome.py | 8.9 |
| scripts/comparative_analysis/09_figures/fig11b_pangenome_plot.py | 4.6 |
| scripts/comparative_analysis/09_figures/fig1_overview.py | 5.1 |
| scripts/comparative_analysis/09_figures/fig1_overview_revised.py | 6.0 |
| scripts/comparative_analysis/09_figures/fig2_chromosome.py | 6.4 |
| scripts/comparative_analysis/09_figures/fig2_chromosome_revised.py | 3.9 |
| scripts/comparative_analysis/09_figures/fig3_shared_plasmid.py | 8.6 |
| scripts/comparative_analysis/09_figures/fig3_shared_plasmid_revised.py | 4.8 |
| scripts/comparative_analysis/09_figures/fig4_bacteriocin_module.py | 7.3 |
| scripts/comparative_analysis/09_figures/fig4_bacteriocin_module_revised.py | 0.0 |
| scripts/comparative_analysis/09_figures/fig5_plasmid2.py | 7.7 |
| scripts/comparative_analysis/09_figures/fig6_phylogeny.py | 10.4 |
| scripts/comparative_analysis/09_figures/fig7_functional.py | 4.7 |
| scripts/comparative_analysis/09_figures/fig8_cazy.py | 4.3 |
| scripts/comparative_analysis/09_figures/fig9_safety.py | 5.4 |
| scripts/comparative_analysis/09_figures/figGA_graphical_abstract.py | 4.4 |
| scripts/comparative_analysis/09_figures/figS1_mapping.py | 6.3 |
| scripts/comparative_analysis/09_figures/figS2_nisin.py | 3.1 |
| scripts/comparative_analysis/09_figures/figS3_candidate_alignment.py | 4.0 |
| scripts/comparative_analysis/09_figures/figS4_pl2_annotation.py | 4.0 |
| scripts/comparative_analysis/09_figures/fig_qc.py | 2.8 |
| scripts/comparative_analysis/09_figures/fig_style.py | 3.5 |
| scripts/comparative_analysis/10_supplement/build_supplementary.py | 5.5 |
| scripts/comparative_analysis/11_manuscript/build_docx.py | 4.0 |
| scripts/comparative_analysis/11_manuscript/build_docx_final.py | 11.5 |
| scripts/comparative_analysis/11_manuscript/build_docx_v2.py | 5.5 |
| scripts/comparative_analysis/12_functional_annotation/parse_functional.py | 5.5 |
| scripts/comparative_analysis/run_all.ps1 | 4.3 |
| scripts/longread_validation/00_manifest/scripts/read_audit2.py | 4.6 |
| scripts/longread_validation/01_references/scripts/build_references.py | 14.1 |
| scripts/longread_validation/01_references/scripts/verify_structure.py | 11.1 |
| scripts/longread_validation/02_mapping/scripts/bam_convert_depth.py | 4.7 |
| scripts/longread_validation/02_mapping/scripts/map_reads2.py | 4.3 |
| scripts/longread_validation/03_supporting_reads/scripts/junction_support2.py | 12.4 |
| scripts/longread_validation/03_supporting_reads/scripts/per_read_models2.py | 5.5 |
| scripts/longread_validation/03_supporting_reads/scripts/per_read_models3.py | 5.9 |
| scripts/longread_validation/04_controls/scripts/copy_number_corrected.py | 4.5 |
| scripts/longread_validation/04_controls/scripts/depth_controls3.py | 9.1 |
| scripts/longread_validation/04_controls/scripts/state_and_junction_counts.py | 4.3 |
| scripts/longread_validation/04_controls/scripts/unit_1bp_test.py | 7.0 |
| scripts/longread_validation/04_controls/scripts/unit_1bp_test2.py | 7.6 |
| scripts/longread_validation/04_controls/scripts/unit_1bp_test3.py | 6.9 |
| scripts/longread_validation/05_figures/scripts/fig_s3.py | 13.4 |
| scripts/longread_validation/06_tables/scripts/build_tables2.py | 10.2 |
| scripts/longread_validation/07_report/scripts/update_manuscript2.py | 9.6 |
| software_versions.txt | 0.9 |

## 5b. 已核实的非敏感命中（扫描假阳性，记录备查）

- `results/tables/derived_longread_validation/junction_candidates2.tsv` 中出现数字串 24879：为比对该表的**数据值**（数值列），非路径/用户名（该表不含任何路径列，已经逐列核对）。
- 功能注释表中 "secretion"（分泌）等词在任何机密类关键词的字面匹配下会命中 "secret" 前缀，属假阳性；已在扫描中用语边界排除。
- 二进制图件（PDF/PNG/SVG）中随机字节可能偶然拼出敏感词；本次已用 PDF 文本层与元数据逐项复核（仅 Matplotlib 生成信息，无个人或 AI 署名）。
