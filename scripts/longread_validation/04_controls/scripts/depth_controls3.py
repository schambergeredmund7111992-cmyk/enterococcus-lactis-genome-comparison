# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
depth_controls3.py — 深度统计（numpy 向量化 bootstrap 版，替代过慢的逐循环实现）
算法1: samtools depth（primary-only，-Q 0 / -Q 20 两档，读自 depth_samtools_regions.tsv）
算法2: 本脚本 CIGAR 解析 —— 同时产出 primary-only 与 all-alignments 两套逐位深度
输出：
  04_controls/depth_region_stats3.tsv      区域统计（算法2，primary 与 all 两套）
  04_controls/depth_ratio_bootstrap3.tsv   相对 chr08 的深度比 + read-resampling bootstrap 95% CI
  04_controls/multimapping_stats.tsv       多映射统计
  04_controls/mixed_structure_counts.tsv   结构群体计数
"""
import gzip
import os
import re
from collections import defaultdict, Counter

import numpy as np

ATTJ = r'<PROJECT_ROOT>\longread_validation'
MAP = os.path.join(ATTJ, '02_mapping')
SR = os.path.join(ATTJ, '03_supporting_reads')
OUT = os.path.join(ATTJ, '04_controls')
LOG = open(os.path.join(ATTJ, '00_manifest', 'logs', 'depth_controls3.log'), 'a', encoding='utf-8')

CIG_RE = re.compile(r'(\d+)([MIDNSHP=X])')
L = {'chr06_full': 2676906, 'chr08_full': 2638187, 'pl2_asassembled': 77437, 'pl1_smbu08': 128837}
S_pos = 1255782; B = 38719; A = 146
j_att_mid = S_pos + A // 2

db_prim = {k: np.zeros(v + 1, dtype=np.int32) for k, v in L.items()}
db_all = {k: np.zeros(v + 1, dtype=np.int32) for k, v in L.items()}
spans_prim = {k: [] for k in L}   # (r0,r1,read, AS)
spans_all = {k: [] for k in L}
n_align = n_prim = n_sec = n_sup = 0
mapq_hist = Counter()
read_refs = defaultdict(set)

with gzip.open(os.path.join(MAP, 'full_db.sam.gz'), 'rt', encoding='utf-8', errors='replace') as f:
    for line in f:
        if line[0] == '@':
            continue
        p = line.rstrip('\n').split('\t')
        ref = p[2]
        if ref == '*':
            continue
        flag = int(p[1]); mapq = int(p[4])
        is_sec = bool(flag & 0x100); is_sup = bool(flag & 0x800)
        n_align += 1
        if is_sec: n_sec += 1
        if is_sup: n_sup += 1
        prim = not is_sec and not is_sup
        if prim:
            n_prim += 1
            mapq_hist[mapq // 10 * 10] += 1
            read_refs[p[0]].add(ref)
        if ref not in L:
            continue
        ops = [(int(a), b) for a, b in CIG_RE.findall(p[5])]
        r0 = int(p[3]) - 1
        pos = r0
        f_m = l_m = None
        tgt = db_prim[ref] if prim else None
        for n, o in ops:
            if o in 'M=X':
                if f_m is None: f_m = pos
                l_m = pos + n
                db_all[ref][pos] += 1; db_all[ref][pos + n] -= 1
                if prim:
                    tgt[pos] += 1; tgt[pos + n] -= 1
                pos += n
            elif o in 'DN':
                pos += n
        if f_m is not None:
            AS = 0
            for t in p[11:]:
                if t.startswith('AS:i:'):
                    AS = int(t[5:])
            spans_all[ref].append((f_m, l_m, p[0]))
            if prim:
                spans_prim[ref].append((f_m, l_m, p[0]))
dep_prim = {k: np.cumsum(db_prim[k][:-1]) for k in L}
dep_all = {k: np.cumsum(db_all[k][:-1]) for k in L}
LOG.write('alignments=%d primary=%d secondary=%d supp=%d\n' % (n_align, n_prim, n_sec, n_sup))

negpos = {}
for line in open(os.path.join(ATTJ, '01_references', 'reference_construction.tsv'), encoding='utf-8').read().splitlines()[1:]:
    p = line.split('\t')
    if p[0].startswith('neg_ctrl_chr08'):
        negpos[p[0]] = int(re.search(r'pos (\d+)', p[3]).group(1))

regions = [
    ('chr08_whole', 'chr08_full', [(0, L['chr08_full'])]),
    ('chr06_whole', 'chr06_full', [(0, L['chr06_full'])]),
    ('pl2_whole', 'pl2_asassembled', [(0, L['pl2_asassembled'])]),
    ('pl1_whole', 'pl1_smbu08', [(0, L['pl1_smbu08'])]),
    ('attJ_±2kb(chr08)', 'chr08_full', [(j_att_mid - 2000, j_att_mid + 2000)]),
    ('attR_±2kb(chr06)', 'chr06_full', [(S_pos + B - 2000, S_pos + B + 2000)]),
    ('attL_±2kb(chr06)', 'chr06_full', [(S_pos - 2000, S_pos + 2000)]),
    ('retained_X_internal(chr06)', 'chr06_full', [(S_pos + A, S_pos + B)]),
    ('dimer_forward_±2kb(pl2)', 'pl2_asassembled', [(B - 2000, B + 2000)]),
    ('dimer_back_wrap_±2kb(pl2)', 'pl2_asassembled', [(L['pl2_asassembled'] - 2000, L['pl2_asassembled']), (0, 2000)]),
    ('monomer_mid_±2kb(pl2)', 'pl2_asassembled', [(B // 2 - 2000, B // 2 + 2000)]),
    ('leftflank_5kb(chr06)', 'chr06_full', [(S_pos - 5000, S_pos)]),
]
for nm, pos in negpos.items():
    if pos:
        regions.append(('%s_±2kb(chr08)' % nm, 'chr08_full', [(pos - 2000, pos + 2000)]))

rows = []
for name, ref, ivs in regions:
    d1 = np.concatenate([dep_prim[ref][a:b] for a, b in ivs])
    d2 = np.concatenate([dep_all[ref][a:b] for a, b in ivs])
    q1 = np.percentile(d1, [25, 50, 75]); q2 = np.percentile(d2, [25, 50, 75])
    rows.append([name, ref, len(d1), round(float(d1.mean()), 1), int(q1[1]), '%d-%d' % (q1[0], q1[2]),
                 round(float(d2.mean()), 1), int(q2[1]), '%d-%d' % (q2[0], q2[2])])
with open(os.path.join(OUT, 'depth_region_stats3.tsv'), 'w', encoding='utf-8') as f:
    f.write('region\ttarget\tn_positions\tdepth_primary_mean\tdepth_primary_median\tdepth_primary_IQR\t'
            'depth_allaln_mean\tdepth_allaln_median\tdepth_allaln_IQR\n')
    for r_ in rows:
        f.write('\t'.join(str(x) for x in r_) + '\n')

# ---------- bootstrap（numpy 向量化，read 重采样） ----------
def boot_ratio(ref, ivs, n_boot=2000, seed=20261002, use_prim=True):
    src = spans_prim if use_prim else spans_all
    # 每 read 的 (贡献给区域的碱基数, 贡献给 chr08 的碱基数)
    byread_reg = defaultdict(int)
    byread_chr = defaultdict(int)
    reg_len = sum(b - a for a, b in ivs)
    for a, b, rd in src[ref]:
        for ra, rb in ivs:
            lo, hi = max(a, ra), min(b, rb)
            if hi > lo:
                byread_reg[rd] += hi - lo
    for a, b, rd in src['chr08_full']:
        byread_chr[rd] += b - a
    reads = sorted(set(byread_reg) | set(byread_chr))
    idx = {rd: i for i, rd in enumerate(reads)}
    reg_arr = np.zeros(len(reads)); chr_arr = np.zeros(len(reads))
    for rd, v in byread_reg.items():
        reg_arr[idx[rd]] = v
    for rd, v in byread_chr.items():
        chr_arr[idx[rd]] = v
    rng = np.random.default_rng(seed)
    n = len(reads)
    pick = rng.integers(0, n, size=(n_boot, n))
    tot_reg = reg_arr[pick].sum(axis=1) / reg_len
    tot_chr = chr_arr[pick].sum(axis=1) / L['chr08_full']
    ratio = tot_reg / np.maximum(tot_chr, 1e-9)
    lo, hi = np.percentile(ratio, [2.5, 97.5])
    return float(np.mean(ratio)), float(lo), float(hi)

brows = []
for name, ref, ivs in regions:
    if name == 'chr08_whole':
        continue
    d1 = np.concatenate([dep_prim[ref][a:b] for a, b in ivs]).mean()
    d2 = np.concatenate([dep_all[ref][a:b] for a, b in ivs]).mean()
    r1 = float(d1 / dep_prim['chr08_full'].mean())
    r2 = float(d2 / dep_all['chr08_full'].mean())
    m, lo, hi = boot_ratio(ref, ivs, use_prim=True)
    brows.append([name, ref, round(r1, 3), round(r2, 3), round(lo, 3), round(hi, 3)])
with open(os.path.join(OUT, 'depth_ratio_bootstrap3.tsv'), 'w', encoding='utf-8') as f:
    f.write('region\ttarget\tratio_primary_algA\tratio_allaln_algB\tbootstrap_CI95_low\tbootstrap_CI95_high\n')
    for r_ in brows:
        f.write('\t'.join(str(x) for x in r_) + '\n')
LOG.write('bootstrap done\n')

# ---------- 多映射统计 ----------
multi = sum(1 for r, s in read_refs.items() if len(s) > 1)
with open(os.path.join(OUT, 'multimapping_stats.tsv'), 'w', encoding='utf-8') as f:
    f.write('item\tvalue\n')
    for k, v in [('full_db_alignments', n_align), ('full_db_primary', n_prim), ('full_db_secondary', n_sec),
                 ('full_db_supplementary', n_sup), ('reads_total', len(read_refs)),
                 ('reads_primary_mapq0', mapq_hist.get(0, 0)),
                 ('reads_primary_mapq>=20', sum(v for k_, v in mapq_hist.items() if k_ >= 20)),
                 ('reads_with_primary_to_>1_ref', multi),
                 ('mapq_hist_primary', dict(sorted(mapq_hist.items())))]:
        f.write('%s\t%s\n' % (k, v))

# ---------- 结构群体 ----------
cls_path = os.path.join(SR, 'read_classification3.tsv')
jr = defaultdict(Counter)
if os.path.exists(cls_path):
    for line in open(cls_path, encoding='utf-8').read().splitlines()[1:]:
        p = line.split('\t')
        for jt in p[1].split(','):
            jr[jt][p[11]] += 1
with open(os.path.join(OUT, 'mixed_structure_counts.tsv'), 'w', encoding='utf-8') as f:
    f.write('# 候选 reads 按 junction 类型 × 模型判定（read_classification3.tsv）\n')
    f.write('junction_type\tstrong\tambiguous\tweak\n')
    for jt, cc in sorted(jr.items()):
        f.write('%s\t%d\t%d\t%d\n' % (jt, cc.get('strong', 0), cc.get('ambiguous', 0), cc.get('weak', 0)))
print('depth_controls3 done')
for r_ in brows:
    print(r_)
