# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
candidate_screen.py — 模块 05：共享质粒细菌素候选模块严谨挖掘
1) 提取两株全部 CDS 翻译
2) BLASTP 全部蛋白 vs 细菌素参考集（E≤1e-3）→ 全命中表
3) 候选位点（GL002641 enterocin P-like / 相邻免疫候选 / GL002679 hiracin JM79-like）：
   - 与最近已知细菌素的对齐级 identity/coverage
   - 前导肽（双甘氨酸 GG）与成熟肽基序（YGNGV 等）扫描
   - ±20 kb 邻域基因列表 + 邻近 MGE 距离
4) pyrodigal 短 ORF 复查（30–150 aa）
输出：search/pl1_vs_bacteriocin_refset.tsv, candidate_table.tsv, neighborhood_*.tsv,
      mge_distance_table.tsv, small_orf_screen.tsv
"""
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
from p3lib import (WORK, GBK06, GBK08, BLASTP, MAKEDB, read_fasta, write_fasta,
                   revcomp, run, open_log)

OUT = os.path.join(WORK, '05_bacteriocins')
SEA = os.path.join(OUT, 'search')
REF = os.path.join(OUT, 'refs')
os.makedirs(SEA, exist_ok=True)
LOG = open_log(os.path.join(OUT, 'logs', 'candidate_screen.log'))

# ---------- 1) GBK → 蛋白 ----------
def gbk_proteins(path):
    """返回 [(locus_tag, seqname, start, end, strand, product, protein_seq)]"""
    feats = []
    cur = None
    qual = None
    seqname = '?'
    in_origin = False
    with open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('LOCUS'):
                seqname = line.split()[1]
                in_origin = False
                continue
            if line.startswith('ORIGIN'):
                if cur:
                    feats.append(cur)
                cur, qual, in_origin = None, None, True
                continue
            if line.startswith('//'):
                in_origin = False
                continue
            if in_origin:
                continue
            m = re.match(r'^     (\S+)\s+(\S+)', line)
            if m:
                # 任意特征起始行：先关闭上一个特征；仅 CDS 才开始收集限定符
                if cur:
                    feats.append(cur)
                if m.group(1) == 'CDS':
                    cur = {'loc': m.group(2), 'q': {}, 'seqname': seqname}
                else:
                    cur = None
                qual = None
                continue
            if cur is not None and line.startswith('                     /'):
                mm = re.match(r'^                     /(\w+)=?(.*)$', line.rstrip('\n'))
                if mm:
                    qual = mm.group(1)
                    cur['q'][qual] = mm.group(2).strip().strip('"')
                continue
            if cur is not None and line.startswith('                     ') and qual:
                cur['q'][qual] += line.strip().strip('"')
    if cur:
        feats.append(cur)
    out = []
    for ft in feats:
        loc = ft['loc']
        strand = '-' if loc.startswith('complement') else '+'
        nums = re.findall(r'(\d+)', loc)
        tr = ft['q'].get('translation', '').rstrip('*')
        if len(nums) >= 2 and tr:
            out.append((ft['q'].get('locus_tag', ''), ft['seqname'], int(nums[0]), int(nums[1]),
                        strand, ft['q'].get('product', ft['q'].get('note', '')), tr))
    return out


prot06 = gbk_proteins(GBK06)
prot08 = gbk_proteins(GBK08)
LOG.write('proteins: smbu06=%d smbu08=%d\n' % (len(prot06), len(prot08)))
write_fasta(os.path.join(SEA, 'smbu06_proteins.faa'), [('%s %s' % (lt, p[:50]), s) for lt, sq, a, b, st, p, s in prot06])
write_fasta(os.path.join(SEA, 'smbu08_proteins.faa'), [('%s %s' % (lt, p[:50]), s) for lt, sq, a, b, st, p, s in prot08])
P06 = {lt: dict(seq=s, chrom=sq, start=a, end=b, strand=st, product=p) for lt, sq, a, b, st, p, s in prot06}
P08 = {lt: dict(seq=s, chrom=sq, start=a, end=b, strand=st, product=p) for lt, sq, a, b, st, p, s in prot08}

# ---------- 2) BLASTP vs 细菌素参考集 ----------
refset = os.path.join(REF, 'bacteriocin_refset.faa')
db = os.path.join(SEA, 'bacteriocin_db')
run([MAKEDB, '-in', refset, '-dbtype', 'prot', '-out', db], log=LOG, check=True)
r = run([BLASTP, '-query', os.path.join(SEA, 'smbu06_proteins.faa'), '-db', db,
         '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen',
         '-evalue', '1e-3', '-max_target_seqs', '10', '-max_hsps', '5', '-seg', 'no', '-num_threads', '8'],
        log=LOG, check=True)
rows = []
for line in r.stdout.strip().split('\n'):
    if not line:
        continue
    p = line.split('\t')
    qcov = (int(p[7]) - int(p[6]) + 1) / float(p[12]) * 100
    rows.append(p[:2] + p[2:12] + [round(qcov, 1)])
with open(os.path.join(SEA, 'pl1_vs_bacteriocin_refset.tsv'), 'w', encoding='utf-8') as f:
    f.write('qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqcov_pct\n')
    for row in rows:
        f.write('\t'.join(map(str, row)) + '\n')
LOG.write('blastp hits: %d\n' % len(rows))
# 至 pl1 的强命中（E<=1e-5 且 qcov>=70 或 length>=40）
strong = [row for row in rows if float(row[10]) <= 1e-5 and (row[12] >= 70 or int(row[3]) >= 40)]
LOG.write('strong hits (E<=1e-5, qcov>=70%% or alen>=40): %d\n' % len(strong))
for row in strong:
    q = row[0].split()[0]
    info = P06.get(q, {})
    LOG.write('  %s (%s) vs %s: pid=%s alen=%s qcov=%s E=%s\n' % (
        q, info.get('chrom', '?'), row[1], row[2], row[3], row[12], row[10]))

# ---------- 3) 候选位点分析 ----------
CANDIDATES = {
    'GL002641': ('smbu06GL002641', 'enterocin P-like (类IIa)'),
    'GL002679': ('smbu06GL002679', 'hiracin JM79-like (类IIa)'),
}
# 相邻候选免疫蛋白：641 上下游 200 bp 内的基因
c641 = P06.get('smbu06GL002641')
neigh_641 = [lt for lt, d in P06.items()
             if d['chrom'] == c641['chrom'] and abs(min(d['start'], c641['start']) - max(d['end'], c641['end'])) < 600
             and lt != 'smbu06GL002641']
LOG.write('near GL002641 (<=600bp): %s\n' % [(lt, P06[lt]['product'][:40]) for lt in neigh_641])
for lt in neigh_641:
    strongly = [row for row in strong if row[0].split()[0] == lt]
    CANDIDATES[lt] = ('%s (adjacent to GL002641)' % lt, 'immunity candidate' if not strongly else 'bacteriocin homolog')

cand_rows = []
for key, (lt, note) in CANDIDATES.items():
    d = P06.get(lt)
    if not d:
        LOG.write('missing %s\n' % lt)
        continue
    aa = d['seq']
    hits = [row for row in rows if row[0].split()[0] == lt]
    best = hits[0] if hits else None
    # 双甘氨酸（GG）前导肽基序
    gg = [m.start() + 1 for m in re.finditer('GG', aa[:50])]
    ygngv = aa.find('YGNGV') + 1 if 'YGNGV' in aa else (aa.find('YDNGI') + 1 if 'YDNGI' in aa else 0)
    # 保守类IIa基序扫描（YGNGV / YDNGI 的宽松变体 YxNGx）
    yxngx = re.search('Y[AVLIMF][NS]G[VA]', aa)
    cand_rows.append(dict(
        locus_tag=lt, chrom=d['chrom'], start=d['start'], end=d['end'], strand=d['strand'],
        length_aa=len(aa), product=d['product'], note=note,
        best_subject=(best[1] if best else ''), best_pident=(best[2] if best else ''),
        best_alen=(best[3] if best else ''), best_evalue=(best[10] if best else ''),
        best_qcov=(best[12] if best else ''), gg_positions=';'.join(map(str, gg)),
        ygngv_pos=ygngv, yxngx_motif=(yxngx.group(0) + '@%d' % (yxngx.start() + 1) if yxngx else '')))
    LOG.write('candidate %s: %daa; best=%s pid=%s alen=%s qcov=%s E=%s; GG=%s; YGNGV@%d; %s\n' % (
        lt, len(aa), cand_rows[-1]['best_subject'][:50], cand_rows[-1]['best_pident'],
        cand_rows[-1]['best_alen'], cand_rows[-1]['best_qcov'], cand_rows[-1]['best_evalue'],
        cand_rows[-1]['gg_positions'], ygngv, cand_rows[-1]['yxngx_motif']))

with open(os.path.join(SEA, 'candidate_table.tsv'), 'w', encoding='utf-8') as f:
    cols = list(cand_rows[0].keys())
    f.write('\t'.join(cols) + '\n')
    for r_ in cand_rows:
        f.write('\t'.join(str(r_[c]) for c in cols) + '\n')

# 邻域（±20 kb）
def neighborhood(lt, tag, window=20000):
    d = P06[lt]
    lo, hi = d['start'] - window, d['end'] + window
    out = []
    for lt2, d2 in P06.items():
        if d2['chrom'] == d['chrom'] and d2['start'] >= lo and d2['end'] <= hi:
            out.append((d2['start'], lt2, d2['strand'], len(d2['seq']), d2['product']))
    out.sort()
    path = os.path.join(SEA, 'neighborhood_%s.tsv' % tag)
    with open(path, 'w', encoding='utf-8') as f:
        f.write('start\tlocus_tag\tstrand\tlength_aa\tproduct\tis_candidate\n')
        for st, lt2, strd, ln, prod in out:
            f.write('%d\t%s\t%s\t%d\t%s\t%s\n' % (st, lt2, strd, ln, prod, 'yes' if lt2 == lt else ''))
    LOG.write('neighborhood %s: %d genes in +/-%dkb\n' % (tag, len(out), window // 1000))
    return path


for tag, (lt, _) in CANDIDATES.items():
    if lt in P06:
        neighborhood(lt, tag)

# MGE 距离表
mge_file = os.path.join(WORK, '04_plasmidome', 'annotation', 'mge_inventory.tsv')
mge_rows = []
with open(mge_file, encoding='utf-8') as f:
    next(f)
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) >= 9 and p[1] == 'Plasmid1' and p[0] == 'smbu06':
            mge_rows.append(dict(lt=p[2], start=int(p[3]), end=int(p[4]), product=p[6], nr=p[7]))
with open(os.path.join(SEA, 'mge_distance_table.tsv'), 'w', encoding='utf-8') as f:
    f.write('candidate\tmge_locus\tmge_start\tmge_end\tdistance_bp\tmge_description\n')
    for tag, (lt, _) in CANDIDATES.items():
        if lt not in P06:
            continue
        c = P06[lt]
        for m in mge_rows:
            if m['end'] < c['start']:
                dist = c['start'] - m['end']
            elif m['start'] > c['end']:
                dist = m['start'] - c['end']
            else:
                dist = 0
            f.write('%s\t%s\t%d\t%d\t%d\t%s\n' % (lt, m['lt'], m['start'], m['end'], dist, (m['nr'] or m['product'])[:80]))
    # 模块区 MGE 密度
    if 'smbu06GL002641' in P06 and 'smbu06GL002679' in P06:
        lo = P06['smbu06GL002641']['start'] - 5000
        hi = P06['smbu06GL002679']['end'] + 5000
        n_in = sum(1 for m in mge_rows if m['start'] >= lo and m['end'] <= hi)
        span = (hi - lo) / 1000.0
        f.write('#module_region\t%d-%d\t%d\tMGE_count=%d\tdensity=%.2f/kb\n' % (lo, hi, int(span * 1000), n_in, n_in / span))
        pl1_len = 128837
        f.write('#plasmid_whole\t1-%d\t%d\tMGE_count=%d\tdensity=%.2f/kb\n' % (
            pl1_len, pl1_len, len(mge_rows), len(mge_rows) / (pl1_len / 1000.0)))
        LOG.write('module region MGE density: %d in %.1f kb = %.2f/kb; whole plasmid: %d in %.1f kb = %.2f/kb\n' % (
            n_in, span, n_in / span, len(mge_rows), pl1_len / 1000.0, len(mge_rows) / (pl1_len / 1000.0)))

# ---------- 4) 短 ORF 复查（pyrodigal） ----------
import pyrodigal
pl1 = read_fasta(WORK + r'\03_replicon_compare\alignments\smbu06_pl1.fa')[0][1]
gf = pyrodigal.GeneFinder(meta=False, min_gene=90)
gf.train(pl1)
genes = gf.find_genes(pl1)
ann_cds = [(d['start'], d['end']) for lt, d in P06.items() if d['chrom'] == 'Plasmid1']
small = []
for i, g in enumerate(genes):
    aa = g.translate()
    if 30 <= len(aa) <= 150:
        covered = any(not (g.end < a or g.begin > b) for a, b in ann_cds)
        gg = 'GG' in aa[:45]
        yxngx = bool(re.search('Y[AVLIMF][NS]G[VA]', aa))
        small.append((g.begin, g.end, g.strand, len(aa), covered, gg, yxngx, aa[:60]))
n_notcov = sum(1 for s in small if not s[4])
LOG.write('pyrodigal pl1: ORFs=%d; small(30-150aa)=%d; not_in_company_annotation=%d; GG-motif=%d\n' % (
    len(genes), len(small), n_notcov, sum(1 for s in small if s[5])))
write_fasta(os.path.join(SEA, 'small_orf_proteins.faa'),
            [('pl1_orf_%d_%d_%s_%daa%s' % (s[0], s[1], s[2], s[3], '_GG' if s[5] else ''), s[7])
             for s in small])
with open(os.path.join(SEA, 'small_orf_screen.tsv'), 'w', encoding='utf-8') as f:
    f.write('start\tend\tstrand\tlength_aa\tin_company_annotation\tdouble_glycine_motif\tclassIIa_motif\tseq60\n')
    for s in small:
        f.write('%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\n' % (s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7]))
print('candidate_screen done')
