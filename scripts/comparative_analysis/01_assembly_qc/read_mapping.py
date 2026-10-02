# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
read_mapping.py — 模块 01：原始 reads 回贴核验（覆盖度/广度/混样）
工具：minimap2 2.31（-x sr 用于 Illumina；-x map-ont 用于 Nanopore；bwa/mem2 无 Windows 构建，
      以 minimap2 sr preset 等效替代并在方法中说明）
参考库：组合库 = smbu06 全组装 + smbu08 全组装 + L. lactis 14B4 chr + MG1363 chr
        → 同株判定 + 乳酸菌污染监测（smbu06 与 194 相同，无需重复纳入 194）
输出：
  01_assembly_qc/mapping/mapping_stats.tsv        每样本每平台：reads、映射率、补充比对
  01_assembly_qc/mapping/coverage_per_replicon.tsv 每复制子平均深度/广度/低覆盖区
  01_assembly_qc/mapping/contamination_summary.tsv 命中对象分类
  01_assembly_qc/mapping/depth_<sample>_<platform>.npz 每碱基深度数组（供作图）
"""
import gzip
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
import numpy as np
from p3lib import (WORK, IN06, IN08, MM2, read_fasta, write_fasta, open_log)

OUT = os.path.join(WORK, '01_assembly_qc', 'mapping')
os.makedirs(OUT, exist_ok=True)
LOG = open_log(os.path.join(WORK, '01_assembly_qc', 'logs', 'read_mapping.log'))

SAMPLES = {
    'smbu06': {
        'illumina': [IN06 + r'\smbu06\1.Cleandata\smbu06.IS350_Clean.1.fq.gz',
                     IN06 + r'\smbu06\1.Cleandata\smbu06.IS350_Clean.2.fq.gz'],
        'nanopore': IN06 + r'\smbu06\1.Cleandata\smbu06.filtered_reads.fq.gz',
        'asm06': True,
    },
    'smbu08': {
        'illumina': [IN08 + r'\smbu08\1.Cleandata\smbu08.IS350_Clean.1.fq.gz',
                     IN08 + r'\smbu08\1.Cleandata\smbu08.IS350_Clean.2.fq.gz'],
        'nanopore': IN08 + r'\smbu08\1.Cleandata\smbu08.filtered_reads.fq.gz',
        'asm06': False,
    },
}
REF14B4 = WORK + r'\02_taxonomy\refs\GCF_003176835.1\ncbi_dataset\data\GCF_003176835.1\GCF_003176835.1_ASM317683v1_genomic.fna'
REFMG = WORK + r'\02_taxonomy\refs\GCF_000009425.1\ncbi_dataset\data\GCF_000009425.1\GCF_000009425.1_ASM942v1_genomic.fna'

# ---------- 组合参考库 ----------
repl = {}   # name -> (seq, sample_tag)
for strain, asm in [('smbu06', IN06 + r'\smbu06\2.Assembly\smbu06.Complete.genome.fasta'),
                    ('smbu08', IN08 + r'\smbu08\2.Assembly\smbu08.Complete.genome.fasta')]:
    for n, s in read_fasta(asm):
        tag = strain
        name = '%s|%s' % (tag, n.split()[0])
        repl[name] = s
for tag, path in [('L_lactis_14B4', REF14B4), ('L_lactis_MG1363', REFMG)]:
    for n, s in read_fasta(path):
        repl['%s|%s' % (tag, n.split()[0])] = s

combo = os.path.join(OUT, 'combined_ref.fa')
write_fasta(combo, list(repl.items()))
LOG.write('combined ref: %d sequences, %d bp total\n' % (len(repl), sum(len(s) for s in repl.values())))
offset = {}
pos = 0
order = []
for k, s in repl.items():
    offset[k] = pos
    pos += len(s)
    order.append(k)
ends = {k: (offset[k], offset[k] + len(repl[k])) for k in order}
total_len = pos


def classify(target):
    if target.startswith('smbu06|'):
        return 'smbu06'
    if target.startswith('smbu08|'):
        return 'smbu08'
    if target.startswith('L_lactis'):
        return 'Lactococcus_ref'
    return target


def parse_cigar(cigar):
    ops = []
    num = ''
    for ch in cigar:
        if ch.isdigit():
            num += ch
            continue
        ops.append((int(num), ch))
        num = ''
    return ops


def process_paf(paf, sample, platform):
    """返回 (rows: dict, depth: np.array[total_len])"""
    n_reads = set()
    n_primary = 0
    n_supp = 0
    by_target = {}
    diff = np.zeros(total_len + 1, dtype=np.int32)
    aln_reads = set()
    with open(paf) as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            qname, qlen, qs, qe, strand, tname = p[0], int(p[1]), int(p[2]), int(p[3]), p[4], p[5]
            ts, te = int(p[7]), int(p[8])
            n_reads.add(qname)
            tp = 'P'
            for field in p[12:]:
                if field.startswith('tp:A:'):
                    tp = field[5]
                    break
            if tp == 'S':
                n_supp += 1
                continue
            n_primary += 1
            aln_reads.add(qname)
            tgt = classify(tname)
            by_target[tgt] = by_target.get(tgt, 0) + 1
            # CIGAR
            cigar = None
            for field in p[12:]:
                if field.startswith('cg:Z:'):
                    cigar = field[5:]
                    break
            if cigar is None:
                # 无CIGAR时按块覆盖
                a, b = offset[tname] + ts, offset[tname] + te
                diff[a] += 1
                diff[b] -= 1
                continue
            ref = offset[tname] + ts
            for ln, op in parse_cigar(cigar):
                if op in 'M=X':
                    diff[ref] += 1
                    diff[ref + ln] -= 1
                    ref += ln
                elif op == 'D':
                    ref += ln          # 删除：参考位不计覆盖（保守）
                elif op == 'N':
                    ref += ln
                # I/S/H 不消耗参考
    depth = np.cumsum(diff[:-1])
    rows = dict(sample=sample, platform=platform, reads_total=len(n_reads),
                reads_primary=n_primary, reads_supp=n_supp,
                reads_aligned=len(aln_reads),
                mapped_pct=round(len(aln_reads) / max(1, len(n_reads)) * 100, 3),
                by_target=';'.join('%s:%d' % (k, v) for k, v in sorted(by_target.items())))
    return rows, depth


stats_rows = []
cov_rows = []
cont_rows = []
for sample, cfg in SAMPLES.items():
    for platform in ['illumina', 'nanopore']:
        paf = os.path.join(OUT, '%s_%s.paf' % (sample, platform))
        if platform == 'illumina':
            cmd = [MM2, '-x', 'sr', '-c', '--secondary=no', '-t', '8',
                   combo, cfg['illumina'][0], cfg['illumina'][1]]
        else:
            cmd = [MM2, '-x', 'map-ont', '-c', '--secondary=no', '-t', '8',
                   combo, cfg['nanopore']]
        LOG.write('>>> %s %s\n' % (sample, platform))
        with open(paf, 'w') as f:
            r = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True,
                               encoding='utf-8', errors='replace')
        LOG.write('<<< rc=%d %s\n' % (r.returncode, (r.stderr or '')[-500:]))
        rows, depth = process_paf(paf, sample, platform)
        stats_rows.append(rows)
        LOG.write('stats: %s\n' % rows)
        # 每复制子深度
        for name in order:
            a, b = ends[name]
            d = depth[a:b]
            if len(d) == 0:
                continue
            mean = d.mean()
            if name.startswith(sample + '|') or True:
                cov_rows.append(dict(sample=sample, platform=platform, target=name,
                                     length=len(d), mean_depth=round(float(mean), 1),
                                     breadth_ge1x=round(float((d >= 1).mean() * 100), 2),
                                     breadth_ge10x=round(float((d >= 10).mean() * 100), 2),
                                     frac_below_quarter_meancov=round(float((d < mean / 4).mean() * 100), 3)))
        np.savez_compressed(os.path.join(OUT, 'depth_%s_%s.npz' % (sample, platform)),
                            depth=depth, total_len=total_len, order=np.array(order))
        # 污染分类（只关心 Lactococcus 命中与未映射）
        cont_rows.append(dict(sample=sample, platform=platform,
                              own=rows['by_target'],
                              lactococcus_hits=sum(v for k, v in [(x.split(':')[0], int(x.split(':')[1]))
                                                                  for x in rows['by_target'].split(';') if x]
                                                   if k == 'Lactococcus_ref'),
                              reads_unmapped=rows['reads_total'] - rows['reads_aligned']))
        print(sample, platform, 'done:', rows['mapped_pct'], '%')

with open(os.path.join(OUT, 'mapping_stats.tsv'), 'w', encoding='utf-8') as f:
    f.write('sample\tplatform\treads_total\treads_primary\treads_supplementary\treads_aligned\tmapped_pct\tby_target\n')
    for r in stats_rows:
        f.write('\t'.join(str(r[k]) for k in ['sample', 'platform', 'reads_total', 'reads_primary',
                                              'reads_supp', 'reads_aligned', 'mapped_pct', 'by_target']) + '\n')
with open(os.path.join(OUT, 'coverage_per_replicon.tsv'), 'w', encoding='utf-8') as f:
    f.write('sample\tplatform\ttarget\tlength\tmean_depth\tbreadth_ge1x\tbreadth_ge10x\tfrac_below_quarter_meancov\n')
    for r in cov_rows:
        f.write('\t'.join(str(r[k]) for k in ['sample', 'platform', 'target', 'length', 'mean_depth',
                                              'breadth_ge1x', 'breadth_ge10x', 'frac_below_quarter_meancov']) + '\n')
with open(os.path.join(OUT, 'contamination_summary.tsv'), 'w', encoding='utf-8') as f:
    f.write('sample\tplatform\tlactococcus_ref_hits\treads_unmapped\ttarget_breakdown\n')
    for r in cont_rows:
        f.write('%s\t%s\t%s\t%s\t%s\n' % (r['sample'], r['platform'], r['lactococcus_hits'],
                                          r['reads_unmapped'], r['own']))
LOG.write('DONE\n')
print('read_mapping done')
