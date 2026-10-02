# FINAL_STATUS — attJ 与串联二聚体长读段验证

**日期**: 2026-10-02 ｜ 输入未被修改 ｜ 工作目录: `<PROJECT_ROOT>\longread_validation\`（ASCII 路径别名 `<PROJECT_ROOT>\longread_validation`）

---

## 结论（顶部一句话）

> **Nanopore long reads are consistent with the proposed attJ and circular tandem-dimer structure, but they additionally reveal that the sequenced smbu08 population is a mixture of two chromosome states in which the retained (element-present) form predominates (~85%) over the excised attJ state (~15%); the delivered smbu08 chromosome assembly represents the minority (excised) state. The circular dimer versus monomer circle could not be independently resolved, and the inter-unit 1-bp difference was not validated.**

（对应规格书证据等级 **B**：attJ 与质粒 junction 均有跨越支持；但单体型/二聚体/环状构型不可判别。**且新增一项规格书未预期的重大发现：样本为混合群体**。）

---

## 1. 各关键 junction 的严格支持 read 数

判据（分析前固定）：单条 read 的**一条 primary 比对连续跨越** junction、mapQ ≥ 20、两侧匹配锚定 ≥ 3 kb（宽松为 ≥ 1 kb）。

| junction | 严格 reads | 宽松 reads | 备注 |
|---|---|---|---|
| attJ（chr08，切离态） | **26–27** | 46 | 其中 41 条在竞争模型分析中 strong（可判别保留态） |
| 质粒串联边界（unit1→unit2） | **992** | 1,883 | 环状分子两条拷贝边界序列相同；旋转帧内 2,910 条跨越比对 |
| 质粒回接边界（origin） | 同上（与正向边界同序列） | 同上 | unit2→unit1 与 unit1→unit2 是同一旋转点的两个拷贝边界 |
| 保留态接合（chr06 flank\|attL） | 149 | 257 | 样本多数状态（规格书未预期） |
| 阴性对照（GC 匹配随机位点 ×3） | 193 / 142 / 248 | — | 高深度普通位置的基线 |
| 打乱对照 | 0 | 0 | |

## 2. 最强 read（attJ；按锚定排序）

- `765ab0bd-139d-48a5-a2fa-eb9ba1365e4d`：长度 21,511 bp；锚定 4,984 | 4,985 bp；mapQ 60；identity 99.61%
- `c3ad006f-f202-45ba-b4ec-3bda61199d48`：长度 14,974 bp；锚定 4,987 | 4,979 bp；mapQ 60；identity 99.72%
- 质粒边界最强 read 例：`df8c6f9c-21c1-4189-bbf8-6b86f810d266`（76,125 bp）、`a7e25edc-5aaf-4400-8cd1-d49f8948fb42`（52,906 bp）

## 3. 竞争模型判定（单体 / 线性二聚体 / 环状二聚体）

1,678 条候选 reads × 10 模型逐一比对（预注册：strong = ΔAS ≥ 10、identity ≥ 90%、覆盖 ≥ 30%）：

- **attJ 态：41 strong / 7 ambiguous**——最佳模型 chr08_full / attJ_junction，次优 chr06_full（保留态），ΔAS 3,632–28,336 → **保留态替代被明确排除**。
- **质粒两条边界：0 strong / 218+202 ambiguous**——全部与单体环/保留态**等分（ΔAS = 0）**。原因：两拷贝除 1 bp 外完全相同，边界序列与单体环任意点、保留态 X 区**序列同构**，跨越 reads 层面**不可能判别**（这不是数据不足，而是结构性不可识别）。
- **1 bp 拷贝间差异（唯一序列判别点）**：read 内配对检验仅 2–4 条 read 可用，δ ∈ {−1, 0, +2}，**不足以验证**。
- 因此：**环状串联二聚体 vs 单体环（双拷贝数）无法由当前数据判别**；线性二聚体（=组装原样线性表示）与环状二聚体的差异仅在回接位点，同样不可判别。
- 组装假象排除：边界序列由 1,883 条（严格 992）独立长 reads 跨越支持（含 2,910 条旋转帧跨越比对），**非 junction 级拼接假象**；高深度（Plasmid2 5.67× 染色体，CI 5.51–5.83）与修正双计后的元件丰度（4.0×，CI 3.94–4.05）支持高拷贝丰度，但不构成自主复制证据。

## 4. 可以现在写进正文的英文句子

> Long-read mappings support the attJ junction state (41 model-distinguishing reads; the retained-state alternative is rejected with ΔAS 3.6–28.3 kb) and place 1,883 reads (≥1 kb anchors, mapQ ≥ 20; 992 strict) across the plasmid's tandem boundary, confirming that the assembled 77.4-kb component is not a junction-level assembly artefact. The tandem-dimer and monomer-circle interpretations are sequence-equivalent except for a single 1-bp inter-unit indel that the available reads cannot resolve.

> Long reads show that the sequenced smbu08 population contains both chromosome states, with the retained (element-present) form predominating (257 vs 46 junction-spanning reads, ≈5.6:1); the delivered smbu08 chromosome assembly therefore represents the minority excised haplotype. The retained-state reads carry the smbu08-defining SNP allele, excluding cross-contamination from smbu06/194.

## 5. 不能写进正文的句子与原因

| 不可写 | 原因 |
|---|---|
| "the circular tandem-dimer state was **directly validated / confirmed** by junction-spanning reads" | 单体型/二聚体序列同构，跨越 reads 不可判别（本报告 §6/§3） |
| "the element is **absent from the smbu08 chromosome**" | 只对交付组装成立；样本 ~85% 为保留态（§9） |
| "**autonomously replicating** / high-copy plasmid **replicates autonomously**" | 高覆盖不证明复制机制（本报告 §7） |
| "**induced** / phage-plasmid" | 无诱导或噬菌体颗粒证据 |
| "两种单元间差异（1 bp）**已由 reads 证实**" | read 内检验统计力不足（n=2–4） |

## 6. 交付文件（绝对路径）

**主报告**
- `<PROJECT_ROOT>\longread_validation\07_report\FINAL_STATUS.md`（本文件）
- `…\07_report\long_read_validation_report.md`（中文）／`…\long_read_validation_report_en.md`（英文）

**图（PNG 600dpi / PDF / SVG）**
- `…\05_figures\png\FigS3_junction_reads.png`、`…\05_figures\pdf\...pdf`、`…\05_figures\svg\...svg`

**表**
- `…\06_tables\Supplementary_Table_S3_junction_read_support.xlsx` 与 `.tsv`
- `…\06_tables\junction_reference_manifest.tsv`
- `…\06_tables\mapping_model_comparison.tsv`
- `…\06_tables\ambiguous_and_negative_reads.tsv`
- `…\06_tables\chromosome_state_and_copy_number.tsv`
- `…\04_controls\state_junction_counts.tsv`、`depth_region_stats3.tsv`、`depth_ratio_bootstrap3.tsv`、`depth_samtools_regions.tsv`、`copy_number_corrected.tsv`、`unit_1bp_withinread*.tsv`、`multimapping_stats.tsv`、`mixed_structure_counts.tsv`

**映射与证据**
- `…\02_mapping\`（full/panel 的 SAM.gz、PAF.gz、sorted BAM + bai、flagstat、idxstats、bam_subsets/±2kb 已索引）
- `…\03_supporting_reads\`（junction_candidates2.tsv、junction_support_summary.tsv、alignment_evidence_*.tsv、supporting_reads_*.fasta、read_classification3.tsv、model_scores3.tsv、paf_subsets/*.paf.gz）
- 被取代的先前会话中间文件移至 `…\_superseded_prior_session\`（含截断的旧映射输出）
- `…\01_references\`（references/*.fasta、reference_construction.tsv、independent_verification.tsv、att_structure_evidence.tsv）
- `…\00_manifest\`（input_inventory.tsv、nanopore_read_stats.tsv、tool_versions.txt、logs/、scripts/、tools/）

**论文更新稿**
- `<WORKDIR>\manuscript_fin_修订投稿版_长读段验证后.docx`（另存，不覆盖原稿）与同名 `.md`
