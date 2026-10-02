# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
parse_functional.py — 模块 12（借鉴参考论文方法）：功能注释与安全分析解析
来源：华大交付注释表（两株各一套）
  KEGG   General_Gene_Annotation/<s>.kegg.list.KEGG2Gene.xls（Level1/Level2 计数 + 基因表）
  COG    General_Gene_Annotation/<s>.cog.list.class.catalog.xls（一级/二级分类 + 基因数）
  GO     General_Gene_Annotation/<s>.go.xls（Ontology/Class/计数）
  CAZy   Pathogen_Analysis/Plant/<s>.cazy.statis_5class.stat.xls + cazy.statis_allclass.stat.xls + cazy.list.anno.xls
  CARD   Pathogen_Analysis/Animal/<s>.card.rgi.anno.xls（Strict/Pass 命中）
  ARDB   Pathogen_Analysis/Animal/<s>.ardb.list.anno.xls（含 Min_Identity 阈值列）
  VFDB   Pathogen_Analysis/Animal/<s>.vfdb.list.anno.xls
  T3SS   Pathogen_Analysis/<s>.effectiveT3.std.stat.xls
输出：12_functional_annotation/summary/*.tsv
"""
import os
import sys
import csv

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
from p3lib import WORK, IN06, IN08, open_log

OUT = os.path.join(WORK, '12_functional_annotation')
SUM = os.path.join(OUT, 'summary')
os.makedirs(SUM, exist_ok=True)
LOG = open_log(os.path.join(OUT, 'logs', 'parse_functional.log'))

ROOTS = {'smbu06': IN06 + r'\smbu06', 'smbu08': IN08 + r'\smbu08'}


def rd(path):
    rows = []
    with open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.rstrip('\n').rstrip('\r')
            if line.strip():
                rows.append(line.split('\t'))
    return rows


def w(name, header, rows):
    with open(os.path.join(SUM, name), 'w', encoding='utf-8') as f:
        f.write('\t'.join(header) + '\n')
        for r in rows:
            f.write('\t'.join(str(x) for x in r) + '\n')
    LOG.write('%s: %d rows\n' % (name, len(rows)))


kegg1, kegg2, cog, go, cazy5, cazyfam, card, ardb, vfdb, t3 = ([] for _ in range(10))
for s, root in ROOTS.items():
    G = root + r'\4.Genome_Function\General_Gene_Annotation'
    A = root + r'\4.Genome_Function\Pathogen_Analysis'
    # KEGG
    for r in rd(G + r'\%s.kegg.list.KEGG2Gene.xls' % s)[1:]:
        if len(r) >= 4:
            kegg1.append([s, r[0], sum(int(x[2]) for x in [r]) if False else '', ''])
            kegg2.append([s, r[0], r[1], int(r[2])])
    # COG
    for r in rd(G + r'\%s.cog.list.class.catalog.xls' % s)[1:]:
        if len(r) >= 5:
            cog.append([s, r[0], r[1], r[2], int(r[3])])
    # GO
    for r in rd(G + r'\%s.go.xls' % s)[1:]:
        if len(r) >= 3:
            go.append([s, r[0], r[1], int(r[2])])
    # CAZy 5 类
    for r in rd(A + r'\Plant\%s.cazy.statis_5class.stat.xls' % s)[1:]:
        if len(r) >= 7 and r[0].strip():
            cazy5.append([s, int(r[1]), int(r[2]), int(r[3]), int(r[4]), int(r[5]), int(r[6])])
    # CAZy 家族（allclass）
    for r in rd(A + r'\Plant\%s.cazy.statis_allclass.stat.xls' % s)[1:]:
        if len(r) >= 2 and r[0].strip():
            try:
                cazyfam.append([s, r[0].strip(), int(r[1])])
            except ValueError:
                pass
    # CARD
    for r in rd(A + r'\Animal\%s.card.rgi.anno.xls' % s)[1:]:
        if len(r) >= 10:
            card.append([s] + r[:10])
    # ARDB
    for r in rd(A + r'\Animal\%s.ardb.list.anno.xls' % s)[1:]:
        if len(r) >= 8:
            ardb.append([s] + r[:8])
    # VFDB
    for r in rd(A + r'\Animal\%s.vfdb.list.anno.xls' % s)[1:]:
        if len(r) >= 8:
            vfdb.append([s, r[0], r[1], r[2], r[3], r[4], r[5], r[7] if len(r) > 7 else ''])
    # T3SS
    for r in rd(A + r'\%s.effectiveT3.std.stat.xls' % s)[1:]:
        if len(r) >= 4:
            t3.append([s, r[1].replace(',', ''), r[2], r[3]])

# KEGG Level1 汇总
agg = {}
for s, l1, l2, n in kegg2:
    agg[(s, l1)] = agg.get((s, l1), 0) + n
w('kegg_level1.tsv', ['strain', 'level1', 'n_genes'], [[s, l1, n] for (s, l1), n in sorted(agg.items())])
w('kegg_level2.tsv', ['strain', 'level1', 'level2', 'n_genes'], kegg2)
w('cog.tsv', ['strain', 'first_class', 'cog_code', 'description', 'n_genes'], cog)
w('go_level2.tsv', ['strain', 'ontology', 'go_class', 'n_genes'], go)
w('cazy_5class.tsv', ['strain', 'AA', 'CBM', 'CE', 'GH', 'GT', 'PL'], cazy5)
w('cazy_families.tsv', ['strain', 'family', 'n_genes'], cazyfam)
w('card_hits.tsv', ['strain', 'gene', 'cut_off', 'pass_bitscore', 'best_bitscore', 'best_aro',
                    'best_identities', 'aro', 'drug_class', 'resistance_mechanism'], card)
w('ardb_hits.tsv', ['strain', 'gene', 'identity', 'evalue', 'subject', 'resistance_type',
                    'min_identity', 'resistance_requirement'], ardb)
w('vfdb_hits.tsv', ['strain', 'gene', 'identity', 'evalue', 'vf_id', 'gi', 'type', 'description'], vfdb)
w('t3ss_summary.tsv', ['strain', 'total_proteins', 'predicted_secreted', 'not_secreted'], t3)

# VFDB 分类计数
vagg = {}
for r in vfdb:
    cat = r[7].split(']')[0].lstrip('[').strip() if r[7] else 'unknown'
    vagg[(r[0], cat)] = vagg.get((r[0], cat), 0) + 1
w('vfdb_category_counts.tsv', ['strain', 'vf_category', 'n_hits'],
  [[s, c, n] for (s, c), n in sorted(vagg.items(), key=lambda x: -x[1])])
print('parse_functional done; KEGG L2 rows:', len(kegg2), '| CARD:', len(card), '| ARDB:', len(ardb), '| VFDB:', len(vfdb))
