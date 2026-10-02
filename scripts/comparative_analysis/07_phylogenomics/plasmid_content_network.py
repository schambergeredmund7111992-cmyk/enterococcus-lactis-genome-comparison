# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
plasmid_content_network.py — 模块 07：共享质粒基因内容（Jaccard）距离网络
原因（写入正文/补方法）：面板中无 ≥99% 同一性的完整质粒（最高 96.7%），backbone 基因树不可靠；
采用 gene-content Jaccard 距离 + average linkage 树（NEIGHBOR-JOINING 等价替代，说明理由）。
1) tblastn 共享质粒 145 蛋白 vs 面板组合库 → 每基因组存在的质粒基因集
2) Jaccard 距离矩阵 → scipy linkage → newick（供 Fig6 使用）
输出：07_phylogenomics/trees/plasmid_gene_content_matrix.tsv / plasmid_jaccard_linkage.nwk
"""
import csv
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
import numpy as np
from p3lib import WORK, TBLASTN, read_fasta, write_fasta, run, open_log

OUT = os.path.join(WORK, '07_phylogenomics')
TRD = os.path.join(OUT, 'trees')
os.makedirs(TRD, exist_ok=True)
LOG = open_log(os.path.join(OUT, 'logs', 'plasmid_content.log'))
AN = os.path.join(WORK, '06_public_panel', 'analysis')

# pl1 蛋白
prot = dict((n.split()[0], s) for n, s in read_fasta(os.path.join(WORK, '05_bacteriocins', 'search', 'smbu06_proteins.faa')))
genes = []
with open(os.path.join(WORK, '04_plasmidome', 'annotation', 'pl1_genes.tsv'), encoding='utf-8') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        if r['region'] == 'pl1_smbu06' and r['locus_tag'] in prot:
            genes.append(r['locus_tag'])
qfa = os.path.join(TRD, 'pl1_gene_queries.faa')
write_fasta(qfa, [(g, prot[g]) for g in genes])
LOG.write('pl1 genes as queries: %d\n' % len(genes))

# 面板去重后的 accession → 其 contig 前缀匹配
meta = [m['accession'] for m in csv.DictReader(open(os.path.join(WORK, '06_public_panel', 'metadata', 'panel_metadata.tsv'), encoding='utf-8'), delimiter='\t')]

r = run([TBLASTN, '-query', qfa, '-db', os.path.join(AN, 'panel_db'),
         '-outfmt', '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen',
         '-evalue', '1e-5', '-seg', 'no', '-max_target_seqs', '2000', '-max_hsps', '1', '-num_threads', '8'],
        log=LOG, check=True)
pres = {}
for line in r.stdout.strip().split('\n'):
    if not line:
        continue
    p = line.split('\t')
    acc = p[1].split('|')[0]
    qcov = (int(p[5]) - int(p[4]) + 1) / float(p[10]) * 100   # 12 列输出:qstart=4,qend=5,qlen=10
    if float(p[2]) >= 70 and qcov >= 70:
        pres.setdefault(acc, set()).add(p[0].split('|')[0])
# 加项目两株（全部基因）
pres['smbu06_project'] = set(genes)
pres['smbu08_project'] = set(genes)
LOG.write('genomes with plasmid gene content: %d\n' % len(pres))

accs = sorted(pres)
G = len(accs)
mat = np.zeros((G, G))
for i in range(G):
    for j in range(i + 1, G):
        a, b = pres[accs[i]], pres[accs[j]]
        inter = len(a & b)
        union = len(a | b)
        d = 1 - inter / union if union else 1.0
        mat[i, j] = mat[j, i] = d

with open(os.path.join(TRD, 'plasmid_gene_content_matrix.tsv'), 'w', encoding='utf-8') as f:
    f.write('accession\tn_genes\tdistance_to_smbu06\n')
    i0 = accs.index('smbu06_project')
    for i, a in enumerate(accs):
        f.write('%s\t%d\t%.4f\n' % (a, len(pres[a]), mat[i, i0]))

from scipy.cluster.hierarchy import linkage, to_tree
Z = linkage(mat[np.triu_indices(G, 1)], method='average')


def to_newick(node, parent_dist):
    if node.is_leaf():
        name = accs[node.id]
        return '%s:%.5f' % (name.replace(' ', '_'), parent_dist - node.dist)
    l = to_newick(node.get_left(), node.dist)
    r = to_newick(node.get_right(), node.dist)
    s = '(%s,%s)' % (l, r)
    if node.dist > 0:
        s += ':%.5f' % (parent_dist - node.dist)
    return s


tree = to_tree(Z)
nwk = to_newick(tree, 0.0) + ';'
with open(os.path.join(TRD, 'plasmid_gene_content.nwk'), 'w') as f:
    f.write(nwk + '\n')
LOG.write('plasmid content tree: %d leaves; distances to smbu06: %s\n' % (
    G, sorted('%s:%.3f' % (accs[i][:15], mat[i, i0]) for i in range(G))[:10]))
print('plasmid_content_network done: %d genomes' % G)
