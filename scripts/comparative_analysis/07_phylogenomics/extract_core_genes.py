# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
extract_core_genes.py — 模块 07：核心标记基因提取（面板 + 两株 + 参考）
标记集：本菌株注释中所有 ribosomal protein 基因 + 11 个看家基因（NR 描述检索），约 70 个
方法（沿用前期已验证流程）：面板基因组 + 两株 + 参考 → 单一 BLAST 库（contig 名加前缀）→
      一次 tblastn（E≤1e-10）→ 每基因组每标记最佳命中 → 抽核酸区段、按 qstart 校框、翻译、
      剔除内部终止 → 每标记一个 FAA 文件（供 MAFFT）
输出：07_phylogenomics/markers/marker_<tag>.faa, core_gene_table.tsv
"""
import os
import sys
import glob
import csv

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
from p3lib import (WORK, IN06, IN08, BLAST, TBLASTN, MAKEDB, read_fasta, write_fasta,
                   revcomp, translate_dna, run, open_log)

OUT = os.path.join(WORK, '07_phylogenomics')
MKD = os.path.join(OUT, 'markers')
os.makedirs(MKD, exist_ok=True)
LOG = open_log(os.path.join(OUT, 'logs', 'extract_core_genes.log'))

# ---------- 1. 选标记基因 ----------
ft_path = os.path.join(WORK, '04_plasmidome', 'annotation', 'full_gene_table.tsv')
rows = []
with open(ft_path, encoding='utf-8') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        rows.append(r)
# 用 NR 描述补充（可选）
nr06 = {}
p = IN06 + r'\smbu06\4.Genome_Function\General_Gene_Annotation\smbu06.nr.list.anno.xls'
with open(p, encoding='utf-8', errors='replace') as f:
    next(f)
    for line in f:
        c = line.rstrip('\n').split('\t')
        if len(c) >= 5 and c[0] not in nr06:
            nr06[c[0]] = c[4]

HOUSEKEEP = ['dnaa', 'dnak', 'groel', 'reca', 'gyrase subunit a', 'gyrase subunit b',
             "rna polymerase subunit alpha", "rna polymerase subunit beta",
             "rna polymerase subunit beta'", 'elongation factor tu', 'elongation factor g']
markers = []
for r in rows:
    if r['strain'] != 'smbu06':
        continue
    t = (r['product'] + ' ' + nr06.get(r['locus_tag'], '')).lower()
    if 'ribosomal protein' in t or any(h in t for h in HOUSEKEEP):
        markers.append(r['locus_tag'])
markers = sorted(set(markers))
LOG.write('标记基因: %d 个\n' % len(markers))

# query 蛋白（smbu06 GBK 翻译）
prot = dict((n.split()[0], s) for n, s in read_fasta(os.path.join(WORK, '05_bacteriocins', 'search', 'smbu06_proteins.faa')))
qfa = os.path.join(MKD, '_marker_queries.faa')
write_fasta(qfa, [(m, prot[m]) for m in markers if m in prot])

# ---------- 2. 合并 BLAST 库 ----------
panel_files = []
for fna in glob.glob(os.path.join(WORK, '06_public_panel', 'downloads', '*', '**', '*.fna'), recursive=True):
    parts = fna.replace('/', os.sep).split(os.sep)
    try:
        i = parts.index('downloads')
        acc = parts[i + 1]
    except (ValueError, IndexError):
        acc = os.path.basename(os.path.dirname(fna))
    panel_files.append((acc, fna))
LOG.write('panel genome files: %d\n' % len(panel_files))

combined = os.path.join(OUT, 'panel_all.fna')
n_seq = 0
with open(combined, 'w') as fo:
    for acc, fna in sorted(panel_files):
        for n, s in read_fasta(fna):
            fo.write('>%s|%s\n%s\n' % (acc, n.split()[0], s))
            n_seq += 1
    for tag, fa in [('smbu06_project', IN06 + r'\smbu06\2.Assembly\smbu06.Complete.genome.fasta'),
                    ('smbu08_project', IN08 + r'\smbu08\2.Assembly\smbu08.Complete.genome.fasta'),
                    ('E_lactis_194', WORK + r'\02_taxonomy\refs\GCF_056582645.1\ncbi_dataset\data\GCF_056582645.1\GCF_056582645.1_ASM5658264v1_genomic.fna'),
                    ('L_lactis_14B4', WORK + r'\02_taxonomy\refs\GCF_003176835.1\ncbi_dataset\data\GCF_003176835.1\GCF_003176835.1_ASM317683v1_genomic.fna'),
                    ('L_lactis_MG1363', WORK + r'\02_taxonomy\refs\GCF_000009425.1\ncbi_dataset\data\GCF_000009425.1\GCF_000009425.1_ASM942v1_genomic.fna'),
                    ('L_lactis_IL1403', WORK + r'\02_taxonomy\refs\GCF_000006865.1\ncbi_dataset\data\GCF_000006865.1\GCF_000006865.1_ASM686v1_genomic.fna')]:
        for n, s in read_fasta(fa):
            fo.write('>%s|%s\n%s\n' % (tag, n.split()[0], s))
            n_seq += 1
LOG.write('combined db sequences: %d\n' % n_seq)
db = os.path.join(OUT, 'panel_all_db')
run([MAKEDB, '-in', combined, '-dbtype', 'nucl', '-out', db], log=LOG, check=True)

# ---------- 3. tblastn ----------
r = run([TBLASTN, '-query', qfa, '-db', db,
         '-outfmt', '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen',
         '-evalue', '1e-10', '-seg', 'no', '-max_target_seqs', '2000', '-max_hsps', '1',
         '-num_threads', '8'], log=LOG, check=True)
seqs = None   # 懒加载主体序列


def get_seq(idx_name):
    global seqs
    if seqs is None:
        seqs = {}
        for n, s in read_fasta(combined):
            seqs[n] = s
    return seqs[idx_name]


best = {}   # (marker, genome) -> hit
for line in r.stdout.strip().split('\n'):
    if not line:
        continue
    p = line.split('\t')
    q, s = p[0], p[1]
    genome = s.split('|')[0]
    evalue = float(p[8])
    key = (q, genome)
    if key not in best or evalue < best[key][0]:
        best[key] = (evalue, p)

table_rows = []
per_marker = {}
for (m, genome), (ev, p) in sorted(best.items()):
    qlen, slen = int(p[10]), int(p[11])
    sstart, send = int(p[6]), int(p[7])
    strand = 1 if send > sstart else -1
    s = get_seq(p[1])
    seg = s[sstart - 1:send] if strand == 1 else revcomp(s[send - 1:sstart])
    # 按 qstart 校框
    qstart = int(p[4])
    off = (qstart - 1) % 3
    seg = seg[off:]
    aa = translate_dna(seg, to_stop=True)
    qcov = (int(p[5]) - qstart + 1) / qlen * 100
    ok = ('*' not in aa) and len(aa) >= 0.7 * qlen / 3 and aa
    table_rows.append((m, genome, p[2], round(qcov, 1), ev, len(aa), 'ok' if ok else 'rejected'))
    if ok:
        per_marker.setdefault(m, []).append((genome, aa))

with open(os.path.join(OUT, 'core_gene_table.tsv'), 'w', encoding='utf-8') as f:
    f.write('marker\tgenome\tpident\tqcov_pct\tevalue\taa_len\tstatus\n')
    for row in table_rows:
        f.write('\t'.join(str(x) for x in row) + '\n')

n_full = 0
n_genomes = len(set(g for _, g, *_ in table_rows))
for m, lst in per_marker.items():
    write_fasta(os.path.join(MKD, 'marker_%s.faa' % m), [('%s|%s' % (g, m), aa) for g, aa in sorted(lst)])
    if len(lst) >= 0.8 * 209:
        n_full += 1
LOG.write('生成 %d 个标记比对文件；覆盖率>=80%%基因组 的标记 %d 个\n' % (len(per_marker), n_full))
print('extract_core_genes done: markers=%d files' % len(per_marker))
