# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
annotate_regions.py — 模块 04：前噬菌体区 + Plasmid2 + 共享质粒的功能注释与 MGE 清单
输出：
  prophage_annotation_full.tsv   前噬菌体 51 CDS + NR 命中 + 功能模块分类
  pl2_annotation_full.tsv        Plasmid2 102 CDS + NR 命中 + 模块分类 + 拷贝归属(unit1/unit2)
  phage_module_table.tsv         模块级汇总
  mge_inventory.tsv              两株全基因组 MGE 相关位点(转座酶/整合酶/噬菌体)
  pl1_key_genes.tsv              共享质粒上的复制/转运/其他关键基因（接合元件检索为阴性）
"""
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
from p3lib import WORK, IN06, IN08, open_log

ANN = os.path.join(WORK, '04_plasmidome', 'annotation')
LOG = open_log(os.path.join(WORK, '04_plasmidome', 'logs', 'annotate_regions.log'))

NR06 = IN06 + r'\smbu06\4.Genome_Function\General_Gene_Annotation\smbu06.nr.list.anno.xls'
NR08 = IN08 + r'\smbu08\4.Genome_Function\General_Gene_Annotation\smbu08.nr.list.anno.xls'


def load_nr(path):
    d = {}
    with open(path, encoding='utf-8', errors='replace') as f:
        next(f)
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) >= 5 and p[0] not in d:
                d[p[0]] = dict(identity=p[1], evalue=p[2], subject=p[3], description=p[4])
    return d


NR = {'smbu06': load_nr(NR06), 'smbu08': load_nr(NR08)}
LOG.write('NR table: smbu06=%d smbu08=%d\n' % (len(NR['smbu06']), len(NR['smbu08'])))


def load_tsv(path):
    rows = []
    with open(path, encoding='utf-8') as f:
        cols = f.readline().rstrip('\n').split('\t')
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) == len(cols):
                rows.append(dict(zip(cols, p)))
    return rows


def module_of(product, nrdesc):
    t = (product + ' ' + nrdesc).lower()
    rules = [
        ('lysogeny', ['integrase', 'site-specific', 'excisionase', 'recombinase', 'repressor', 'rghr', 'xis']),
        ('DNA_packaging', ['terminase', 'portal']),
        ('head', ['capsid', 'head', 'prohead', 'scaffold', 'major head']),
        ('tail', ['tail', 'baseplate', 'bppu', 'tape measure', 'structural protein', 'fiber', 'sheath']),
        ('lysis', ['holin', 'amidase', 'endolysin', 'lysin', 'hemolysin', 'xhl', 'autolysin']),
        ('recombination_DNA', ['rect', 'yqa', 'exonuclease', 'recombination', 'nuclease', 'helicase', 'primase', 'polymerase', 'methylase', 'restriction']),
        ('regulation', ['transcriptional', 'regulator', 'sigma', 'antirepressor', 'kil']),
    ]
    for mod, kws in rules:
        if any(k in t for k in kws):
            return mod
    return 'hypothetical'


rows_pro = load_tsv(os.path.join(ANN, 'prophage_genes.tsv'))
rows_pl2 = load_tsv(os.path.join(ANN, 'pl2_genes.tsv'))
rows_pl1 = load_tsv(os.path.join(ANN, 'pl1_genes.tsv'))
LOG.write('loaded: prophage=%d pl2=%d pl1=%d\n' % (len(rows_pro), len(rows_pl2), len(rows_pl1)))

# 附 NR
def attach(rows, strain):
    for r in rows:
        nr = NR[strain].get(r['locus_tag'], {})
        r['nr_subject'] = nr.get('subject', '')
        r['nr_identity'] = nr.get('identity', '')
        r['nr_evalue'] = nr.get('evalue', '')
        r['nr_description'] = nr.get('description', '')
        r['functional_module'] = module_of(r['product'], r['nr_description'])
    return rows


attach(rows_pro, 'smbu06')
attach(rows_pl2, 'smbu08')
attach(rows_pl1, 'smbu06')
# pl1 smbu08 部分
for r in rows_pl1:
    if r['region'] == 'pl1_smbu08':
        nr = NR['smbu08'].get(r['locus_tag'], {})
        r['nr_subject'] = nr.get('subject', '')
        r['nr_identity'] = nr.get('identity', '')
        r['nr_evalue'] = nr.get('evalue', '')
        r['nr_description'] = nr.get('description', '')
        r['functional_module'] = module_of(r['product'], r['nr_description'])

# pl2 拷贝归属
for r in rows_pl2:
    r['copy'] = 'unit1' if int(r['start']) < 38720 else 'unit2'


def write(rows, path):
    cols = ['region', 'chrom', 'copy', 'locus_tag', 'start', 'end', 'strand', 'length_aa',
            'product', 'functional_module', 'nr_subject', 'nr_identity', 'nr_evalue', 'nr_description']
    with open(path, 'w', encoding='utf-8') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')


write(rows_pro, os.path.join(ANN, 'prophage_annotation_full.tsv'))
write([r for r in rows_pl2], os.path.join(ANN, 'pl2_annotation_full.tsv'))

# 模块汇总
from collections import Counter
for tag, rows in [('prophage', rows_pro), ('plasmid2', rows_pl2)]:
    c = Counter(r['functional_module'] for r in rows)
    LOG.write('%s modules: %s\n' % (tag, dict(c)))
with open(os.path.join(ANN, 'phage_module_table.tsv'), 'w', encoding='utf-8') as f:
    f.write('region\tfunctional_module\tn_genes\tgenes\n')
    for tag, rows in [('prophage_smbu06_38.7kb', rows_pro), ('plasmid2_smbu08_77.4kb', rows_pl2)]:
        groups = {}
        for r in rows:
            groups.setdefault(r['functional_module'], []).append(r)
        for mod, rs in groups.items():
            genes = ';'.join('%s(%saa,%s)' % (r['locus_tag'], r['length_aa'],
                                              (r['nr_description'][:40] or r['product'][:40]))
                             for r in sorted(rs, key=lambda x: int(x['start'])))
            f.write('%s\t%s\t%d\t%s\n' % (tag, mod, len(rs), genes))

# pl1 关键基因（复制/稳定/转运/接合阴性证据）
KEYWORDS = ['replication', 'replicase', 'repA', 'relaxase', 'mob', 'tra', 'conjug', 'transfer',
            'toxin', 'antitoxin', 'partition', 'segregation', 'stability', 'killer', 'immunity',
            'single-stranded', 'helicase', 'primase', 'dnaa', 'resolvase', 'integrase', 'transposase']
key_rows = []
for r in rows_pl1:
    t = (r['product'] + ' ' + r.get('nr_description', '')).lower()
    if any(k in t for k in KEYWORDS):
        key_rows.append(r)
write(key_rows, os.path.join(ANN, 'pl1_key_genes.tsv'))
LOG.write('pl1 key genes: %d\n' % len(key_rows))
for r in key_rows:
    LOG.write('  %s %s %s\n' % (r['locus_tag'], r['product'][:60], r['nr_description'][:60]))

# MGE 清单（全基因组检索）
MGE_KW = ['transposase', 'insertion sequence', 'is256', 'integrase', 'recombinase', 'resolvase',
          'phage', 'prophage', 'terminase', 'conjugal', 'relaxase', 'replication initiator']
mge_rows = []
all_genes = load_tsv(os.path.join(ANN, 'full_gene_table.tsv'))
for r in all_genes:
    t = (r['product'] + ' ' + NR[r['strain']].get(r['locus_tag'], {}).get('description', '')).lower()
    hits = [k for k in MGE_KW if k in t]
    if hits:
        nr = NR[r['strain']].get(r['locus_tag'], {})
        mge_rows.append(dict(strain=r['strain'], chrom=r['chrom'], locus_tag=r['locus_tag'],
                             start=r['start'], end=r['end'], strand=r['strand'],
                             product=r['product'], nr_description=nr.get('description', ''),
                             matched=hits))
with open(os.path.join(ANN, 'mge_inventory.tsv'), 'w', encoding='utf-8') as f:
    f.write('strain\tchrom\tlocus_tag\tstart\tend\tstrand\tproduct\tnr_description\tmatched_keywords\n')
    for r in sorted(mge_rows, key=lambda x: (x['strain'], x['chrom'], int(x['start']))):
        f.write('\t'.join(str(r[k]) for k in ['strain', 'chrom', 'locus_tag', 'start', 'end', 'strand',
                                              'product', 'nr_description']) + '\t' + ';'.join(r['matched']) + '\n')
LOG.write('MGE inventory rows: %d\n' % len(mge_rows))
from collections import Counter as C2
LOG.write('MGE by chrom: %s\n' % dict(C2((r['strain'], r['chrom']) for r in mge_rows)))
print('annotate_regions done; MGE rows:', len(mge_rows))
