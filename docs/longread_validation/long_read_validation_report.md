# smbu08 attJ 与串联二聚体结构的长读段验证报告

**日期**：2026-10-02  
**数据**：smbu08 原始 Nanopore clean reads（`smbu08.filtered_reads.fq.gz`，125,722 reads / 1,217.5 Mb / N50 13,827 bp，SHA-256 `aeafb703…`）+ smbu08/smbu06 完整组装；smbu06 Nanopore reads 作为阳性对照。  
**性质**：纯计算验证，未进行任何湿实验。

---

## 0. 一句话结论

> **Nanopore long reads are consistent with the proposed attJ and circular tandem-dimer structure, but they additionally reveal that the sequenced smbu08 population is a mixture of two chromosome states in which the retained (chr06/194-like, element-present) chromosome predominates (~85%) over the excised attJ state (~15%); the delivered smbu08 chromosome assembly represents the minority (excised) state.**

即：attJ 连接状态真实存在且被可区分模型的 reads 强支持（41 条 strong），质粒串联边界被 ~1,900 条长reads 跨越支持；但单体型（monomer）与二聚体在序列上除 1 bp 外完全同构、无法由跨越 reads 判别；并且样本为**保留态占多数的混合群体**，这一点修正了"smbu08 染色体缺失该元件"的表述——缺的是**组装所捕获的单倍型**，不是样本的主体。

---

## 1. 输入审计（`00_manifest/`）

| 项目 | 值 |
|---|---|
| Nanopore reads | 125,722 条 / 1,217,502,251 bp；mean 9,684；median 6,804；N50 13,827；max 174,754；≥7 kb 占 48.7%；≥30 kb 占 3.5% |
| fq.gz SHA-256 | `aeafb703ef844c1ba2134e3711161f5fdc734ff35e0c0e9501bc6bccdd7911b9` |
| smbu08 染色体 | 2,638,187 bp（sha `fb0c3268…`，与交付组装一致） |
| smbu08 质粒 fasta | 含 Plasmid1 128,837 bp + Plasmid2 77,437 bp（sha `9d846394…`） |
| smbu06 染色体 | 2,676,906 bp（sha `58d48697…`） |
| 工具 | minimap2 2.31-r1302（Windows 构建）、NCBI BLAST+ 2.17、samtools 1.24（MSYS2 ucrt64 官方镜像，`00_manifest/tools/bin/`）、Python 3.13（miniconda base） |

原始 reads 在任务开始时已被归档为 `<DATA_ROOT>\smbu08.tar.gz` 并从 Downloads 删除；本次由该归档恢复至原路径并逐文件核对完整性（127 条目全部解出，组装 SHA-256 与前期记录一致）。**原始文件未被修改。**

## 2. 结构复核（独立第二实现，`01_references/independent_verification.tsv`）

字节级、不依赖任何比对器边界放置：

- `chr08 == chr06 删除 chr06:1,255,783–1,294,501（0-based 1255782–1294501，38,719 bp）后再替换 1 个碱基`：
  - 左翼 1,255,782 bp **0 差异**；右翼 1,382,405 bp **1 个 SNP**（chr08:2,404,144 C ↔ chr06:2,442,863 T）。
- `attL`（chr06 1,255,783–1,255,928）、`attJ`（chr08 同坐标）、`attR`（chr06 1,294,502–1,294,647）：
  - **attJ 与 attR 字节完全相同**；attL 与 attJ/attR 在 146 bp 中差 59 个位点。
- 元件 `attL + X + attR` = 38,865 bp（X = 38,573 bp）。
- `unit1`（38,719 bp）与 `revcomp(attL+X)` **旋转等价（offset 130）**——即质粒单元与切离片段在全长上精确等同（非"相似"）；`unit2` = unit1 删除 unit1 第 36,517 位（poly-A 内 1 bp）。
- 全部 16 条已存参考序列与本次重算逐字节一致；阴性对照来源核实无误。
- 共享质粒 Plasmid1 与 smbu06 的同名质粒旋转等价（offset 17,178）。

## 3. Junction 参考面板（`01_references/reference_construction.tsv`）

13 条参考：chr06/chr08 全染色体、pl1、pl2 组装原样、attJ 窗口、保留态窗口、unit1/unit2、两条 10 kb junction 参考、monomer 环模型（两个线性化点）、旋转帧二聚体（回接内部化）、打乱阴性对照、3 条 GC 匹配随机阴性对照。每条记录构造方法、长度与 SHA-256；阳性/阴性对照均经复核。

## 4. 映射（`02_mapping/`，`00_manifest/logs/map_reads2.log`）

minimap2 `-x map-ont`，全库（chr06+chr08+pl1+pl2）与 panel 库各跑**两次**（`-c`→PAF；`-a --MD --secondary=yes`→SAM），SAM 经 samtools 转 sorted BAM + 索引（`full_db.sorted.bam` 1.30 GB / `panel_db.sorted.bam` 1.15 GB），保留未过滤原始输出、命令、版本与 flagstat/idxstats。

- 全库比对 283,979 条：primary 123,168 / secondary 145,097 / supplementary 15,714；99.11% reads 可映射。
- **92.2% 的 primary 比对 mapQ=0**——因为两株染色体除元件与 1 SNP 外完全相同，染色体 reads 在两参考间等分，无法唯一指派。这是本数据集的核心特征，所有基于 mapQ 的判据都必须显式声明其含义（panel 内 mapQ vs 全库 mapQ）。

## 5. 跨 junction 支持 reads（预注册判据）

**判据**（分析前固定）：单条 read 的一条比对连续跨越 junction（M 覆盖）、primary、mapQ ≥ 20、两侧匹配锚定 ≥ 1 kb（严格 ≥ 3 kb）；SAM 与 PAF 双路径独立解析并对照。

| junction | 宽松 reads | 严格 reads | 说明 |
|---|---|---|---|
| **attJ（chr08，切离态）** | 46 | 26–27 | 有可区分模型的强支持 |
| 保留态左接合（chr06 flank\|attL） | 257 | 149 | 见 §9（样本多数状态） |
| 保留态右接合（chr06 X\|attR） | 154 | 87 | |
| 切离态右接合（chr08 attJ\|flank） | 48 | 27 | |
| **质粒串联边界（u1\|u2）** | **1,883** | **992** | 全库 mapQ≥20；pl2 单帧内 2,908–2,910 条跨越比对 |
| 质粒回接边界（origin） | 同一序列（环状分子两个拷贝边界序列相同），由旋转帧 dimer_back_internal 覆盖：2,910 条跨越比对 | | |
| 阴性对照 ×3（GC 匹配随机位点） | 248 / 142 / 193 | 高深度下普通染色体位置的基线跨越数 | |
| 打乱对照 | 0 | 0 | |

最强 attJ read 例（按两侧锚定排序）：`765ab0bd-139d-48a5-a2fa-eb9ba1365e4d`（长度 21,511 bp，锚定 4,984 | 4,985 bp，mapQ 60，identity 99.61%）；`c3ad006f-f202-45ba-b4ec-3bda61199d48`（14,974 bp，锚定 4,987 | 4,979，mapQ 60，identity 99.72%）。完整清单见 `03_supporting_reads/alignment_evidence_*.tsv` 与 Supplementary Table S3。

## 6. 每条支持 read 的竞争模型比对（`03_supporting_reads/read_classification3.tsv`）

对 1,678 条候选（含二聚体 junction 全部 primary 跨越 reads）逐一对 10 个模型单独映射，取最佳 AS，按预注册阈值判定（strong：ΔAS ≥ 10 且 identity ≥ 90% 且覆盖 ≥ 30%；ambiguous：ΔAS < 10 或 identity < 90%）：

| junction 类别 | strong | ambiguous | 主要最佳模型 |
|---|---|---|---|
| attJ（chr08） | **41** | 7 | chr08_full / attJ_junction（次优 chr06_full，ΔAS 3,632–28,336，identity 97–99%） |
| 保留态（chr06） | 150 | 108 | chr06_full / retained_window |
| 质粒正向边界 | **0** | 218 | unit1_to_unit2_junction / pl2（**ΔAS = 0**） |
| 质粒回接边界 | 0 | 202 | 同上 |
| 阴性对照 | 0 | 952 | chr08_full（普通染色体位置） |

**结论与成因（必须如实写入）**：质粒边界跨越 reads 大量存在，但**全部与单体环/保留态模型等分**——因为串联二聚体的两个拷贝除 1 bp 外完全相同，其"拷贝边界"序列与单体环的任意一点、以及保留态染色体 X 区序列**在本层面上同构**。这不是数据不足，而是**该结构在跨越 reads 层面不可判别**；判定依赖：(i) 组件在组装中的独立存在（77,437 bp contig 由 reads 全程支持）；(ii) 拷贝间 1 bp 差异（§8）；(iii) 深度（§7）。

## 7. 深度、拷贝数与替代解释排除（`04_controls/`）

两种算法：samtools depth（primary-only；`-Q 0` 与 `-Q 20` 两档）与自实现 CIGAR 逐位解析（同时产出 primary-only 与 all-alignments 两套）。pl1 上两实现完全吻合（485.3 vs 485.9×），在重复区域差异可完全由 secondary 计数口径解释。

| 区域 | primary 深度 | all-alignments | 相对 chr08（bootstrap CI95） |
|---|---|---|---|
| chr08 全体 | 177.0 | 377.8 | 1.00 |
| chr06 全体 | 186.0 | 408.4 | 1.05（1.04–1.07） |
| **Plasmid2** | **1,002.9** | 2,644.6 | **5.67（5.51–5.83）** |
| attJ ±2 kb | 73.6 | 153.5 | 0.42（0.34–0.50） |
| attR ±2 kb（保留态独有） | 241.4 | 697.9 | 1.36（1.22–1.51） |
| 元件 X 内部 | 825.0 | 2,444.0 | 4.66（4.49–4.83） |
| 质粒正向边界 ±2 kb | 2,291.9 | 2,792.4 | 12.95（12.45–13.50） |
| 质粒回接边界 ±2 kb | 123.4 | 994.4 | 0.70–2.63（线性化切口所致，非真实低覆盖） |
| 阴性对照 ±2 kb | 148.7–269.9 | 336.2–553.6 | 0.84–1.53（≈染色体水平 ✓） |

**拷贝数（修正串联双计）**：以 read 为单位、对同一 read 的多次比对按 1/n 加权：元件碱基 106.5 Mb vs 染色体碱基 907.7 Mb → 元件丰度 ≈ **4.0 倍染色体**（CI 3.94–4.05；按 77.4 kb 二聚体解释 ≈ 4 拷贝/染色体当量，按 38.7 kb 单体解释 ≈ 8 拷贝）。**高覆盖不等于自主复制**——本数据不能证明复制机制。

**深度均匀性**：染色体存在 2.5 倍级 ori/ter 样梯度（95×–259×），因此所有窗口深度比均需局部归一或同位置对照；junction 计数不受影响（同一位置）。

**多映射**：92.2% primary mapQ=0（见 §4）；元件 reads 在 chr06-X/pl2-u1/pl2-u2 三处等分——拷贝数分析已按此修正。

## 8. 单元间 1 bp 差异的 read 内检验（`04_controls/unit_1bp_*.tsv`）

这是单体型与二聚体**唯一的序列判别点**（位于 poly-A：unit1 位置 36,517）。批量 copy1-vs-copy2 分布比较**不可判别**（read 在两拷贝间的指派由比对分数决定，而分数差恰来自该 1 bp 本身——参考指派偏倚，已在报告中说明）。唯一干净的方法是**同一 read 内配对**（同一 read 同时覆盖两拷贝对应窗口）：

- 三种实现穷举（单帧 primary；pl2+旋转双帧；任意参考泛化定位）：合计仅 **4 条唯一 read** 同时覆盖两拷贝窗口；δ ∈ {+2, −1, 0, 0}，各实现内均值 +1.0 / −0.5 / −0.33（bootstrap CI 均跨 0）；对照 A-run（两拷贝序列本应相同）在全部 30+ 配对中除 1 例外均为 0。**n 太小，无法判别**。
- 结论：**当前数据不足以独立验证该 1 bp 差异**，因此也不足以据此排除"单体环双拷贝数"解释。该 1 bp 的最终判定依赖组装图/厂商流程（不可得）。

## 9. 关键新发现：smbu08 样本为两种染色体状态的混合群体（`04_controls/state_junction_counts.tsv`）

**证据链（多重独立）：**

1. **同位置 junction 计数（Nanopore，同一 read 集、同一判据）**：
   - 保留态接合（chr06 flank\|attL）：**257 条 loose / 149 条 strict reads**
   - 切离态接合（chr08 flank\|attJ）：**46 条 loose / 26 条 strict reads**
   - → 保留:切离 ≈ **5.6 : 1（≈ 85% : 15%）**
2. **不依赖比对器的 20-mer 邻接检验**（read 序列内 flank-mer 与 attL/attJ-mer 距离 ≤ 300 bp）：184 条 read 含"flank→attL"邻接（中位距离 100 bp）；35 条含"flank→attJ"；两集合**互斥**（无一 read 同时含两种——物理上不可能）。
3. **阳性对照 smbu06（应为纯保留）**：保留 195 条 vs 切离 8 条（24:1，背景 ≈ 4%）——方法在纯样本上表现正确；smbu08 的切离信号（46）远高于该背景。
4. **独立技术与独立映射（前期 Illumina/Nanopore 数据）**：attR（保留独有）vs attJ（切离独有）窗口深度——smbu08 Illumina 196.6 vs 82.3；smbu08 Nanopore 193.3 vs 66.5；smbu06 对照 attJ 仅 12.8–29.3（基线）。方向一致，确认 smbu08 存在 smbu06 没有的切离亚群。
5. **排除 smbu06/194 污染**：SNP 位点（chr08:2,404,144）等位计数 C 462 : T 5（1.07% ≈ 测序错误水平）。若"保留态分子"来自 smbu06/194 污染（携 T），在 ~85% 比例下应见 ~85% T——实测 1%，**排除污染**；保留态分子携 smbu08 自身等位，属该样本内部群体。
6. **attR/attJ 窗口 primary 深度**（同数据集）：180.6 vs 61.2（≈3:1），方向一致（窗口法受覆盖梯度影响，绝对值仅作辅助）。

**推论**：交付的 smbu08 染色体组装（切离态）捕获的是**少数单倍型**；样本主体为保留态。可能解释（均为假说，需厂商流程信息才能定论）：组装流程在混合位点解析到了切离路径、或重复/质粒解析步骤将 prophage 从染色体上删除。

## 10. 论文可用/不可用语句

**可用（长读段支持，含限定）**：
> Long-read mappings support the attJ junction state (41 model-distinguishing reads; retained-state alternative rejected with ΔAS 3.6–28.3 kb) and place 1,883 reads (≥1 kb anchors, mapQ ≥ 20; 992 strict) across the plasmid's tandem boundary, confirming the assembled 77.4-kb component is not a junction-level assembly artefact. The tandem-dimer and monomer-circle interpretations are sequence-equivalent except for a single 1-bp inter-unit indel that the available reads cannot resolve.

**必须新增（混合群体）**：
> Long reads show that the sequenced smbu08 population contains both chromosome states, with the retained (element-present) form predominating (~257 vs 46 junction-spanning reads, ≈5.6:1); the delivered smbu08 chromosome assembly therefore represents the minority excised haplotype. The retained-state reads carry the smbu08-defining SNP allele, excluding cross-contamination from smbu06/194.

**不可用（无证据）**：`induction`、`autonomous replication`、`phage-plasmid`、任何将高覆盖等同于自主复制的表述；`directly validated circular dimer`（单体/二聚体未判别）；`the element is absent from the smbu08 chromosome`（应改为"absent from the delivered assembly / minority haplotype"）。

## 11. 交付与复现

- 图：`05_figures/{png,pdf,svg}/FigS3_junction_reads.*`（A–F 六面板；实线=read 证据，虚线=推断）
- 表：`06_tables/Supplementary_Table_S3_junction_read_support.{xlsx,tsv}`、`junction_reference_manifest.tsv`、`mapping_model_comparison.tsv`、`ambiguous_and_negative_reads.tsv`、`chromosome_state_and_copy_number.tsv`
- 读段与证据：`03_supporting_reads/`（候选与支持 read 的 FASTA、逐 read 证据表）、`02_mapping/bam_subsets/`（±2 kb BAM 子集）
- 全部脚本位于各模块 `scripts/`，日志位于 `00_manifest/logs/`；每个数字可追溯到命令与输出文件。
- 环境与第三方工具来源：`00_manifest/tool_versions.txt`（samtools 为 MSYS2 官方镜像预编译包，因 win-64 无 bioconda 构建；依赖闭包一并列入 `00_manifest/tools/`）。

**局限**：无 samtools 官方 Windows 版（用 MSYS2 官方镜像替代）；无组装图，无法判定厂商组装流程在混合位点的路径选择；1 bp 差异与单体型/二聚体判别统计力不足；所有结论为计算推断。
