# smbu06 与 smbu08 比较基因组完整分析报告（v2，含功能注释与安全性分析）

**报告日期**：2026-10-02（v2：按课题组既有论文格式重排；新增功能注释/CAZy/安全性/ANI/泛基因组分析，全部图件内嵌）
**工作区**：`<PROJECT_ROOT>\comparative_analysis\`
**证据分级**：A = 直接由本次组装/reads/逐碱基比对支持；B = 比较基因组推断；C = 现有数据不可判定（不进入结论）。
**全部分析仅基于测序数据与公共数据库；无任何实验数据。**

---

## 一、结论（先读）

1. **两株为同一 194 样谱系的近等基因衍生物（A）**。smbu06 染色体与公开的 *Enterococcus lactis* 194（GCF_056582645.1）**字节级完全相同**（0 SNP）；smbu08 与其仅差**一个 38,719 bp 的 PBSX 样元件缺失 + 1 个远端 SNP**。ANIb 两株对 194 均为 100.00%（图 10）。
2. **共享质粒 128,837 bp 两株完全相同（A，两条独立证据）**：规范化 SHA-256 相等；BLASTN 全覆盖 100.000%、0 SNP、0 InDel（图 3）。质粒无接合元件。
3. **染色体差异 = 一个 PBSX 样可切离元件 + 1 个 SNP（A/B）**（图 2）。元件含整合酶/PBSX 末端酶/门户/卷尺蛋白/尾-基板/裂解盒（图 5、图 S4）；公司前噬菌体工具未检出。
4. **smbu08 的"额外质粒"（Plasmid2，77,437 bp）是该元件切离态（A）**：38,719 bp 单元的二聚体（第 2 拷贝缺 1 bp），深度为染色体的 4.1–5.0 倍（图 5、图 S1）。
5. **细菌素盒（A 同源/B 功能注释）**：enterocin P-like（94.4%，基序变体 YDNGI）＋ Sakacin-A 免疫同源物（10 bp 下游）＋ hiracin JM79-like（59.3%）；显著共现；模块区无 MGE 富集（图 4）。
6. **功能基因组景观（A，注释层面）**：KEGG 以代谢类（碳水化合物代谢最高）为主；COG 为典型乳酸菌谱；GO 以催化活性/结合与代谢过程为主（图 7）。
7. **CAZy 酶谱（A，注释层面）**：GH 86 + GT 44 + CBM 20 + CE 10 + PL 3 + AA 2（165 基因），以 GH13/GH1/GH3/GH32 与 GT2/GT4 为主（图 8）——与"乳源、可广泛利用寡糖并组装多糖链"的基因组特征一致。
8. **安全相关注释（A，如实报告）**：CARD 严格命中 **4 个高一致性耐药基因**（均位于染色体）：**AAC(6′)-Ii 99.45%**、msrC 97.15%、eatAv 96.00%、efrA 82.57%（图 9）；ARDB 命中全部低于库阈值（41–46% vs 80%）；VFDB 128 条全为 "Predicted" 且以管家蛋白为主。**注意**：这与"未检出耐药基因"的表述不符，安全性讨论需按精确数值如实陈述。
9. **面板分布（A/B）**：159 株去重面板中无 ≥80% 覆盖的质粒匹配（exact 仅 194 本身）；最近亲属共享 60.8–67.7% 骨架（96.7–98.5% id）；盒基因 10.7%/5.7%/6.9%，显著共现（OR 24.2, p=1.1e-5；免疫基因不单现 p=1.7e-10）（图 6）。
10. **两株均无 nisin 簇（A）**：双层检索对照通过（F44 阳性 ✓、IL1403 阴性 ✓）（图 S2）。
11. **泛基因组（A）**：smbu06/smbu08/194/IDCC 2105/CX 2-6_2 五基因组共享庞大核心，E. faecium 64/3 仅共享泛种核心；支持近等基因判定（图 11）。
12. **不可判定（C）**：传播方向/时间、切离机制、正式菌株登记——需实验或元数据层面解决。

---

## 二、功能注释与安全性（v2 新增）

### 2.1 KEGG / COG / GO（图 7）
- **KEGG**：代谢类占比最高，其中**碳水化合物代谢**为 Level-2 最富集类别，其次为氨基酸代谢与膜转运——典型的乳源乳酸菌特征。
- **COG**：翻译/核糖体、氨基酸与碳水化合物转运代谢占主要份额，含多个糖摄取 PTS/ABC 系统。
- **GO**：分子功能以催化活性、结合为主；生物过程以代谢过程与细胞过程为主。
- smbu08 与 smbu06 各分类差异 ≤2 基因（近等基因预期）。

图 7 功能基因组景观（KEGG/COG/GO，smbu06；smbu08 基本一致）：

[[FIG:09_figures/Fig7_functional/Fig7.png|15.0]]

### 2.2 CAZy（图 8）
- smbu06：GH 86、GT 44、CBM 20、CE 10、PL 3、AA 2（共 165）；smbu08：GH 87（共 166）。
- 主要家族：GH13、GH1、GH3、GH32、GT2、GT4、CE4 等；CBM 辅助底物捕获。
- 仅注释层面结论；未测定多糖产量（限制见第五节）。

图 8 CAZy 酶谱：

[[FIG:09_figures/Fig8_cazy/Fig8.png|15.0]]

### 2.3 安全性（CARD/ARDB/VFDB，图 9，如实报告）
| 库 | 结果 | 说明 |
|---|---|---|
| CARD（严格） | **4 个高一致性命中（两株相同、均在染色体）** | AAC(6′)-Ii 99.45%（氨基糖苷灭活）；msrC 97.15%、eatAv 96.00%（靶保护）；efrA 82.57%（外排） |
| ARDB | 17 条/株，identity 41–46% | 全部低于库 Min_Identity（80%）→ 不计为耐药决定因子 |
| VFDB | 128 条/株，全部 "Predicted" | 最高 ~78%，以管家蛋白（如 EF-Tu）为主 |

**结论表述建议**："基因组携带数个已知耐药基因（含 AAC(6′)-Ii 99.5%）；未测定其表达与功能。" 不可写"无耐药基因"。

图 9 安全性相关注释：

[[FIG:09_figures/Fig9_safety/Fig9.png|15.0]]

---

## 三、比较基因组学（v2 新增）

### 3.1 ANI 热图（图 10）
smbu06/smbu08/194 三者 ANIb=100.00%；对参考论文所用比较株 **CX 2-6_2 = 98.82%**、IDCC 2105 = 98.29%；对外群 E. faecium/E. faecalis 显著更低（见矩阵）。

[[FIG:09_figures/Fig10_ani/Fig10.png|13.5]]

### 3.2 泛基因组（图 11）
六基因组（smbu06、smbu08、194、IDCC 2105、CX 2-6_2、E. faecium 64/3）pyrodigal 统一预测 + BLASTP(≥60%/60%) 聚类：五株肠球菌共享庞大核心；E. faecium 仅共享泛种核心。核心/附属/特异计数见 `pangenome_summary.tsv`。

[[FIG:09_figures/Fig11_pangenome/Fig11.png|15.0]]

---

## 四、核心证据链（与 v1 相同，可复算）

| # | 事实 | 证据文件 | 图 |
|---|---|---|---|
| 1 | 输入完整性（232 文件 SHA-256） | `00_manifest/input_inventory.tsv` | — |
| 2 | reads 无混样（100% 映射；Lactococcus ≤0.028%） | `01_assembly_qc/mapping/*.tsv` | S1 |
| 3 | smbu06=194 字节级；ANIb 100% | `02_taxonomy/*.tsv` | 10 |
| 4 | 共享质粒完全相同（双证据） | `03_replicon_compare/identity_hashes.tsv` 等 | 3 |
| 5 | 38,719 bp 元件 + 1 SNP | `chr_difference_summary.tsv` 等 | 2 |
| 6 | Plasmid2 = 切离单元二聚体（4.1–5.0×深度） | `pl2_structure_summary.tsv` | 5 |
| 7 | 细菌素盒（94.4%/免疫/59.3%） | `05_bacteriocins/search/*` | 4, S3 |
| 8 | nisin 阴性（对照通过） | `nisin_status_matrix.tsv` 等 | S2 |
| 9 | 面板分布与共现统计 | `06_public_panel/analysis/*`、`08_statistics/*` | 6 |
| 10 | KEGG/COG/GO/CAZy/安全 | `12_functional_annotation/summary/*` | 7,8,9 |
| 11 | 泛基因组 | `12_functional_annotation/summary/pangenome_summary.tsv` | 11 |

## 五、限制（不写进结论的内容）

- 全部为序列计算：无表达/加工/分泌/抑菌/多糖产量实验。
- 面板非随机采样，频率为描述统计。
- 泛基因组采用单连接启发式聚类（≥60%/60%），可能合并旁系家族（簇级计数偏保守）。
- 元件切离为单时点观测；传播方向与时间不可判定。
- **安全性结论须按精确数值表述**（4 个 CARD 高一致性命中），并与既有论文中"无耐药基因"的说法存在出入——建议核对该表述依据。


---

## 六、全部图件索引（内嵌；PNG/PDF/SVG 三种格式见 03_图件/）

图 1 两株闭环复制子概览与组装统计：

[[FIG:09_figures/Fig1_overview/Fig1.png|14.5]]

图 2 染色体比较（一个可切离元件 + 1 个 SNP）：

[[FIG:09_figures/Fig2_chromosome/Fig2.png|14.5]]

图 3 共享质粒 128,837 bp 完全保守（两条独立证据）：

[[FIG:09_figures/Fig3_shared_plasmid/Fig3.png|13.0]]

图 4 细菌素盒（enterocin P-like / 免疫 / hiracin JM79-like）：

[[FIG:09_figures/Fig4_bacteriocin_module/Fig4.png|14.5]]

图 5 smbu08 特异 77.4 kb 复制子 = 切离态前噬菌体：

[[FIG:09_figures/Fig5_plasmid2/Fig5.png|14.5]]

图 6 系统发育与公共面板分布：

[[FIG:09_figures/Fig6_phylogeny/Fig6.png|15.0]]

图 S1 reads 回贴与覆盖度；图 S2 nisin 双层检索对照；图 S3 候选肽比对；图 S4 元件完整注释：

[[FIG:09_figures/FigS1_mapping/FigS1.png|14.0]]

[[FIG:09_figures/FigS2_nisin/FigS2.png|14.0]]

[[FIG:09_figures/FigS3_candidate_alignment/FigS3.png|14.0]]

[[FIG:09_figures/FigS4_pl2_annotation/FigS4.png|14.0]]

图形摘要（Graphical abstract）：

[[FIG:09_figures/GraphicalAbstract/GraphicalAbstract.png|15.0]]
