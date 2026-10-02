# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
stats_analysis.py — 模块 08：分布统计、共现与敏感性分析
输入（模块 06 输出）：
  gene_presence_matrix.tsv（641/642/679，阈值 70/70；60/70；80/80）
  gene_hits_all.tsv（原始命中，含 contig/坐标 → 连锁分析）
  pl1_distribution.tsv（共享质粒分布分层）
  nisin_panel_status.tsv；panel_metadata.tsv（去重后）
输出：08_statistics/…tsv（频数、Fisher、物种分层、连锁、阈值敏感性）
"""
import csv
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
from p3lib import WORK, open_log
from scipy.stats import fisher_exact

OUT = os.path.join(WORK, '08_statistics')
os.makedirs(OUT, exist_ok=True)
LOG = open_log(os.path.join(OUT, 'logs', 'stats_analysis.log'))
PAN = os.path.join(WORK, '06_public_panel', 'analysis')


def load_tsv(path):
    with open(path, encoding='utf-8') as f:
        return list(csv.DictReader(f, delimiter='\t'))


meta = {m['accession']: m for m in load_tsv(os.path.join(WORK, '06_public_panel', 'metadata', 'panel_metadata.tsv'))}
pres = load_tsv(os.path.join(PAN, 'gene_presence_matrix.tsv'))
dist = {d['accession']: d for d in load_tsv(os.path.join(PAN, 'pl1_distribution.tsv'))}
nisin = {n['accession']: n for n in load_tsv(os.path.join(PAN, 'nisin_panel_status.tsv'))}
LOG.write('panel: %d genomes; with pl1 hits: %d\n' % (len(meta), len(dist)))

GENES = ['enterocinP_like', 'hiracinJM79_like', 'sakacinA_immunity']
def matrix(col):
    m = {g: {} for g in GENES}
    for row in pres:
        m[row['gene']][row['accession']] = row[col] == 'yes'
    return m

M = matrix('presence_70_70')
N = len(meta)

# 频数（总体 + 物种分层）
species = {}
for a, m in meta.items():
    sp = m['organism']
    species.setdefault(sp, []).append(a)
rows = []
for g in GENES:
    n = sum(1 for a in meta if M[g].get(a))
    rows.append(('all', g, n, N, round(n / N * 100, 1)))
    for sp, accs in sorted(species.items(), key=lambda x: -len(x[1])):
        if len(accs) < 5:
            continue
        n2 = sum(1 for a in accs if M[g].get(a))
        rows.append((sp, g, n2, len(accs), round(n2 / len(accs) * 100, 1)))
with open(os.path.join(OUT, 'gene_frequency.tsv'), 'w', encoding='utf-8') as f:
    f.write('group\tgene\tn_carriers\tn_total\tpct\n')
    for r in rows:
        f.write('\t'.join(map(str, r)) + '\n')
LOG.write('frequency table rows: %d\n' % len(rows))

# Fisher 共现（641 & 679；641 & immunity）
def fisher(ga, gb, col='presence_70_70'):
    Mm = matrix(col)
    a11 = sum(1 for x in meta if Mm[ga].get(x) and Mm[gb].get(x))
    a10 = sum(1 for x in meta if Mm[ga].get(x) and not Mm[gb].get(x))
    a01 = sum(1 for x in meta if not Mm[ga].get(x) and Mm[gb].get(x))
    a00 = len(meta) - a11 - a10 - a01
    orr, p = fisher_exact([[a11, a10], [a01, a00]], alternative='greater')
    # OR 95% CI（对数法，Haldane-Anscombe 校正）
    import math
    aa, bb, cc, dd = a11 + 0.5, a10 + 0.5, a01 + 0.5, a00 + 0.5
    lor = math.log((aa * dd) / (bb * cc))
    se = math.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    ci = (math.exp(lor - 1.96 * se), math.exp(lor + 1.96 * se))
    return a11, a10, a01, a00, orr, p, ci

f_rows = []
for pair in [('enterocinP_like', 'hiracinJM79_like'), ('enterocinP_like', 'sakacinA_immunity'),
             ('hiracinJM79_like', 'sakacinA_immunity')]:
    for col in ['presence_70_70', 'presence_60_70', 'presence_80_80']:
        a11, a10, a01, a00, orr, p, ci = fisher(pair[0], pair[1], col)
        f_rows.append((pair[0], pair[1], col, a11, a10, a01, a00, round(orr, 3), p,
                       '%.3f-%.3f' % ci))
        LOG.write('Fisher %s & %s @%s: both=%d, onlyA=%d, onlyB=%d, neither=%d; OR=%.2f p=%.3g CI=%.2f-%.2f\n'
                  % (pair[0], pair[1], col, a11, a10, a01, a00, orr, p, ci[0], ci[1]))
with open(os.path.join(OUT, 'cooccurrence_fisher.tsv'), 'w', encoding='utf-8') as f:
    f.write('geneA\tgeneB\tthreshold\ta11_both\ta10_onlyA\ta01_onlyB\ta00_neither\todds_ratio\tfisher_p_greater\tOR_CI95\n')
    for r in f_rows:
        f.write('\t'.join(map(str, r)) + '\n')

# 连锁（同 contig、间距）— 用 641 与 679 的原始命中
hits = load_tsv(os.path.join(PAN, 'gene_hits_all.tsv'))
by = {}
for h in hits:
    if h['gene'] in ('enterocinP_like', 'hiracinJM79_like'):
        by.setdefault(h['accession'], {}).setdefault(h['gene'], []).append(h)
link_rows = []
for acc, d in sorted(by.items()):
    if 'enterocinP_like' in d and 'hiracinJM79_like' in d:
        for h1 in d['enterocinP_like']:
            for h2 in d['hiracinJM79_like']:
                if h1['contig'] == h2['contig']:
                    gap = abs(max(int(h1['sstart']), int(h1['send'])) - min(int(h2['sstart']), int(h2['send'])))
                    link_rows.append((acc, h1['contig'], gap))
with open(os.path.join(OUT, 'linkage_641_679.tsv'), 'w', encoding='utf-8') as f:
    f.write('accession\tcontig\tgap_bp\n')
    for r in link_rows:
        f.write('\t'.join(map(str, r)) + '\n')
if link_rows:
    gaps = [r[2] for r in link_rows]
    LOG.write('linkage: %d genomes, gap min=%d median=%d max=%d\n' % (
        len(set(r[0] for r in link_rows)), min(gaps), sorted(gaps)[len(gaps) // 2], max(gaps)))

# 质粒关联：pl1 携带状态 × 基因
pl1_status = {a: d['class'] for a, d in dist.items()}
cross = []
CLASSES = ['full_length_near_identical', 'full_length_divergent', 'partial_backbone',
           'fragment', 'trace', 'no_plasmid_hit']
for g in GENES:
    for cls in CLASSES:
        n_g = sum(1 for a in meta if (pl1_status.get(a, 'no_plasmid_hit') == cls) and M[g].get(a))
        n_t = sum(1 for a in meta if pl1_status.get(a, 'no_plasmid_hit') == cls)
        cross.append((g, cls, n_g, n_t))
LOG.write('cross classes: %s\n' % {cls: sum(1 for a in meta if pl1_status.get(a, 'no_plasmid_hit') == cls) for cls in CLASSES})
with open(os.path.join(OUT, 'plasmid_gene_association.tsv'), 'w', encoding='utf-8') as f:
    f.write('gene\tpl1_class\tn_with_gene\tn_class_total\n')
    for r in cross:
        f.write('\t'.join(map(str, r)) + '\n')

# nisin 汇总
call = {}
for a, n in nisin.items():
    call[n['cluster_call']] = call.get(n['cluster_call'], 0) + 1
LOG.write('nisin panel: %s\n' % call)
print('stats_analysis done; panel N=%d' % N)
