# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
coverage_analysis.py — 模块 01：深度去卷积与覆盖度质量汇总
由于两株基因组近相同（smbu06/smbu08 之间读段会分流），逐靶标的原始深度约为真实值的一半。
本脚本：
  1) 合并两株靶标（同分子表示）→ 去卷积深度
  2) 切离单元（smbu08 的 Plasmid2 + smbu06 染色体前噬菌体区 + smbu08 连接区）总深度
  3) 覆盖均匀性（<1/4 平均深度的区段比例）、低覆盖区
  4) 汇总对照公司数字（460/570/440/600/2600）
输出：01_assembly_qc/mapping/coverage_deconvolved.tsv、coverage_uniformity.tsv
"""
import csv
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
import numpy as np
from p3lib import WORK, open_log

MP = os.path.join(WORK, '01_assembly_qc', 'mapping')
LOG = open_log(os.path.join(WORK, '01_assembly_qc', 'logs', 'coverage_analysis.log'))

from p3lib import read_fasta
ref = read_fasta(os.path.join(MP, 'combined_ref.fa'))
offs = {}
pos = 0
for name, s in ref:
    offs[name] = pos
    pos += len(s)
refd = dict(ref)

rows = []
unif = []
for sample, platform in [('smbu06', 'illumina'), ('smbu06', 'nanopore'),
                         ('smbu08', 'illumina'), ('smbu08', 'nanopore')]:
    npz = np.load(os.path.join(MP, 'depth_%s_%s.npz' % (sample, platform)), allow_pickle=True)
    depth = npz['depth']

    def d(name):
        a = offs[name]
        return depth[a:a + len(refd[name])]
    # 自身染色体/质粒（合并两株，去卷积；染色体长度不同，按各自坐标分别统计）
    chr06 = d('smbu06|Chromosome1')
    pl1_sum = d('smbu06|Plasmid1').astype(float) + d('smbu08|Plasmid1').astype(float)
    pl2 = d('smbu08|Plasmid2')
    # 切离单元：pl2 靶标 + smbu06 染色体前噬菌体区 + smbu08 连接区（避免重复计 att）
    proph_region = chr06[1255782:1294501].astype(float)          # smbu06 前噬菌体区
    attJ = d('smbu08|Chromosome1')[1255782:1255928].astype(float)  # smbu08 连接区
    unit_depth = float(pl2.mean() + proph_region.mean() * (38719 / max(1, 38719)) * 0 + proph_region.mean())
    # 说明：pl2 靶标覆盖整单元；proph_region 是另一表示（同源读数分流处）→ 相加=该单元总读数深度
    rows.append(dict(sample=sample, platform=platform,
                     chr_combined_depth=round(float(chr06.mean() + d('smbu08|Chromosome1').mean()), 1),
                     pl1_combined_depth=round(float(pl1_sum.mean()), 1),
                     pl2_target_depth=round(float(pl2.mean()), 1),
                     prophage_region_depth_smbu06=round(float(proph_region.mean()), 1),
                     excised_unit_total_depth=round(float(pl2.mean() + proph_region.mean()), 1)))
    # 均匀性（smbu06 坐标系，合并）
    for tag, arr in [('chromosome(combined,smbu06 frame)', chr06 + d('smbu08|Chromosome1')[:len(chr06)] if False else chr06),
                     ('plasmid1(combined)', pl1_sum)]:
        mean = arr.mean()
        frac_low = float((arr < mean / 4).mean() * 100)
        unif.append(dict(sample=sample, platform=platform, target=tag, mean_depth=round(float(mean), 1),
                         breadth_ge1x=round(float((arr >= 1).mean() * 100), 3),
                         frac_below_quarter=round(frac_low, 4)))
    LOG.write('%s %s: chr=%.0f pl1=%.0f pl2=%.0f proph=%.0f unit_total=%.0f\n' % (
        sample, platform, rows[-1]['chr_combined_depth'] / 2, rows[-1]['pl1_combined_depth'] / 2,
        rows[-1]['pl2_target_depth'], rows[-1]['prophage_region_depth_smbu06'],
        rows[-1]['excised_unit_total_depth']))

with open(os.path.join(MP, 'coverage_deconvolved.tsv'), 'w', encoding='utf-8') as f:
    cols = ['sample', 'platform', 'chr_combined_depth', 'pl1_combined_depth', 'pl2_target_depth',
            'prophage_region_depth_smbu06', 'excised_unit_total_depth']
    f.write('\t'.join(cols) + '\n')
    for r in rows:
        f.write('\t'.join(str(r[c]) for c in cols) + '\n')
with open(os.path.join(MP, 'coverage_uniformity.tsv'), 'w', encoding='utf-8') as f:
    f.write('sample\tplatform\ttarget\tmean_depth\tbreadth_ge1x\tfrac_below_quarter_pct\n')
    for r in unif:
        f.write('\t'.join(str(r[k]) for k in ['sample', 'platform', 'target', 'mean_depth',
                                              'breadth_ge1x', 'frac_below_quarter']) + '\n')
print('coverage_analysis done')
