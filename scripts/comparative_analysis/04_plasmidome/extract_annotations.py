# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
extract_annotations.py — 模块 04：从交付 GBK 提取关键区段基因注释
区段：
  A) smbu06 前噬菌体区 1,255,783–1,294,501 (38,719 bp)
  B) smbu08 Plasmid2 (77,437 bp)
  C) smbu06 Plasmid1 / smbu08 Plasmid1（共享质粒,145/146 CDS）
输出：
  04_plasmidome/annotation/prophage_genes.tsv
  04_plasmidome/annotation/pl2_genes.tsv
  04_plasmidome/annotation/pl1_genes.tsv
  04_plasmidome/annotation/replicon_gene_summary.tsv
"""
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
from p3lib import WORK, GBK06, GBK08, read_fasta, ASM06, ASM08, gc_content, open_log

OUT = os.path.join(WORK, '04_plasmidome', 'annotation')
os.makedirs(OUT, exist_ok=True)
LOG = open_log(os.path.join(WORK, '04_plasmidome', 'logs', 'extract_annotations.log'))


def parse_gbk(path):
    """GenBank 解析（多记录支持）：返回特征列表 dict，含 seqname"""
    feats = []
    cur = None
    qual = None
    seqname = '?'
    in_origin = False
    with open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('LOCUS'):
                seqname = line.split()[1] if len(line.split()) > 1 else '?'
                in_origin = False
                continue
            if line.startswith('ORIGIN'):
                if cur:
                    feats.append(cur)
                cur = None
                qual = None
                in_origin = True
                continue
            if line.startswith('//'):
                in_origin = False
                continue
            if in_origin:
                continue
            m = re.match(r'^     (\S+)\s+(\S+)', line)
            if m and m.group(1) in ('CDS', 'tRNA', 'rRNA', 'source', 'gene'):
                if cur:
                    feats.append(cur)
                cur = {'type': m.group(1), 'loc': m.group(2), 'quals': {}, 'seqname': seqname}
                qual = None
                continue
            if cur is not None and line.startswith('                     /'):
                mm = re.match(r'^                     /(\w+)=?(.*)$', line.rstrip('\n'))
                if mm:
                    qual = mm.group(1)
                    cur['quals'][qual] = mm.group(2).strip().strip('"')
                continue
            if cur is not None and line.startswith('                     ') and qual:
                cur['quals'][qual] += line.strip().strip('"')
                continue
    if cur:
        feats.append(cur)
    out = []
    for ft in feats:
        loc = ft['loc']
        strand = -1 if loc.startswith('complement') else 1
        nums = re.findall(r'(\d+)', loc)
        if len(nums) >= 2:
            s, e = int(nums[0]), int(nums[1])
            out.append(dict(type=ft['type'], seqname=ft['seqname'], start=s, end=e,
                            strand='-' if strand < 0 else '+', **ft['quals']))
    return out


def region_table(feats, chrom, lo, hi, label):
    rows = []
    for ft in feats:
        if ft['type'] != 'CDS':
            continue
        if ft.get('seqname') != chrom:
            continue
        if ft['start'] >= lo and ft['end'] <= hi:
            tr = ft.get('translation', '')
            rows.append(dict(region=label, chrom=chrom, locus_tag=ft.get('locus_tag', ''),
                             start=ft['start'], end=ft['end'], strand=ft['strand'],
                             length_nt=ft['end'] - ft['start'] + 1,
                             length_aa=len(tr.rstrip('*')),
                             product=ft.get('product', ft.get('note', '')),
                             db_xref=ft.get('db_xref', '')))
    return rows


def write_tsv(path, rows, cols):
    with open(path, 'w', encoding='utf-8') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join(str(r.get(c, '')).replace('\t', ' ') for c in cols) + '\n')


feats06 = parse_gbk(GBK06)
feats08 = parse_gbk(GBK08)
LOG.write('GBK parsed: smbu06 features=%d, smbu08 features=%d\n' % (len(feats06), len(feats08)))

cols = ['region', 'chrom', 'locus_tag', 'start', 'end', 'strand', 'length_nt', 'length_aa', 'product', 'db_xref']
proph = region_table(feats06, 'Chromosome1', 1255783, 1294501 + 200, 'prophage_1,255,783-1,294,501')
write_tsv(os.path.join(OUT, 'prophage_genes.tsv'), proph, cols)
LOG.write('prophage CDS: %d\n' % len(proph))
for r in proph:
    LOG.write('  %s %s-%s(%s) %saa %s\n' % (r['locus_tag'], r['start'], r['end'], r['strand'], r['length_aa'], r['product'][:70]))

pl2 = region_table(feats08, 'Plasmid2', 1, 77437, 'plasmid2')
write_tsv(os.path.join(OUT, 'pl2_genes.tsv'), pl2, cols)
LOG.write('pl2 CDS: %d\n' % len(pl2))
for r in pl2:
    LOG.write('  %s %s-%s(%s) %saa %s\n' % (r['locus_tag'], r['start'], r['end'], r['strand'], r['length_aa'], r['product'][:70]))

pl1_06 = region_table(feats06, 'Plasmid1', 1, 128837, 'pl1_smbu06')
pl1_08 = region_table(feats08, 'Plasmid1', 1, 128837, 'pl1_smbu08')
write_tsv(os.path.join(OUT, 'pl1_genes.tsv'), pl1_06 + pl1_08, cols)
LOG.write('pl1 CDS: smbu06=%d smbu08=%d\n' % (len(pl1_06), len(pl1_08)))

# 汇总（GC、基因密度）
s_chr = read_fasta(ASM06['Chromosome1'])[0][1]
proph_seq = s_chr[1255782:1294501]
pl2_seq = [r for r in read_fasta(ASM08['Plasmid1']) if 'Plasmid2' in r[0]][0][1]
# 全基因组 CDS 表（供其他模块复用）
all06 = [dict(r, strain='smbu06') for r in
         [dict(type=ft['type'], chrom=ft['seqname'], locus_tag=ft.get('locus_tag', ''),
               start=ft['start'], end=ft['end'], strand=ft['strand'],
               length_nt=ft['end'] - ft['start'] + 1,
               length_aa=len(ft.get('translation', '').rstrip('*')),
               product=ft.get('product', ft.get('note', '')), db_xref=ft.get('db_xref', ''))
          for ft in feats06 if ft['type'] == 'CDS']]
all08 = [dict(r, strain='smbu08') for r in
         [dict(type=ft['type'], chrom=ft['seqname'], locus_tag=ft.get('locus_tag', ''),
               start=ft['start'], end=ft['end'], strand=ft['strand'],
               length_nt=ft['end'] - ft['start'] + 1,
               length_aa=len(ft.get('translation', '').rstrip('*')),
               product=ft.get('product', ft.get('note', '')), db_xref=ft.get('db_xref', ''))
          for ft in feats08 if ft['type'] == 'CDS']]
write_tsv(os.path.join(OUT, 'full_gene_table.tsv'), all06 + all08,
          ['strain', 'chrom', 'locus_tag', 'start', 'end', 'strand', 'length_nt', 'length_aa', 'product', 'db_xref'])
LOG.write('full gene table: smbu06=%d smbu08=%d\n' % (len(all06), len(all08)))

summary = [
    dict(region='prophage_smbu06', length=len(proph_seq), gc_pct=round(gc_content(proph_seq), 2),
         n_cds=len(proph), cds_per_kb=round(len(proph) / len(proph_seq) * 1000, 2),
         coding_density_pct=round(sum(r['length_nt'] for r in proph) / len(proph_seq) * 100, 1),
         chromosome_gc=round(gc_content(s_chr), 2)),
    dict(region='plasmid2_smbu08', length=len(pl2_seq), gc_pct=round(gc_content(pl2_seq), 2),
         n_cds=len(pl2), cds_per_kb=round(len(pl2) / len(pl2_seq) * 1000, 2),
         coding_density_pct=round(sum(r['length_nt'] for r in pl2) / len(pl2_seq) * 100, 1),
         chromosome_gc=round(gc_content(read_fasta(ASM08['Chromosome1'])[0][1]), 2)),
]
write_tsv(os.path.join(OUT, 'replicon_gene_summary.tsv'), summary,
          ['region', 'length', 'gc_pct', 'n_cds', 'cds_per_kb', 'coding_density_pct', 'chromosome_gc'])
for s in summary:
    LOG.write('summary: %s\n' % s)
print('extract_annotations done')
