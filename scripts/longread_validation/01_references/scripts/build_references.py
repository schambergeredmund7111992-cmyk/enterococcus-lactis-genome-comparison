# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
build_references.py — 01_references：从原始 FASTA 独立重建结构并构建 junction 参考面板
（不引用此前项目结论；所有坐标、方向、长度、重复边界均从序列比对自动推导并输出证据表）
输出：
  01_references/references/*.fasta        全部参考序列
  01_references/reference_construction.tsv 每个参考的构造方法/来源/坐标/长度/SHA-256
  01_references/att_structure_evidence.tsv attL/attR/attJ 长度、identity、差异位点、字节比较
  01_references/units_and_junctions.tsv    pl2 两单元边界、junction 坐标、1bp 差异位置
"""
import os
import sys
import hashlib

ROOT = r'<PROJECT_ROOT>\longread_validation'
IN08 = r'<PROJECT_ROOT>\smbu08'
IN06 = r'<PROJECT_ROOT>\smbu06'
MM2 = r'<TOOLS_ROOT>\minimap2-win\minimap2-2.31-r1302-windows-x86_64-ucrt64\minimap2.exe'
BLAST = r'<TOOLS_ROOT>\ncbi-blast-2.17.0+\bin'
sys.path.insert(0, r'<PROJECT_ROOT>\comparative_analysis\00_manifest\scripts')
from p3lib import read_fasta, write_fasta, revcomp, min_rotation, run, sha256_file

REFDIR = os.path.join(ROOT, '01_references', 'references')
os.makedirs(REFDIR, exist_ok=True)
LOG = open(os.path.join(ROOT, '00_manifest', 'logs', 'build_references.log'), 'a', encoding='utf-8')

# ---------- 载入（复制子身份：按长度与组装统计判断，不按文件名猜测） ----------
chr06 = read_fasta(IN06 + r'\smbu06\2.Assembly\smbu06.Chromosome.fasta')[0][1]
chr08 = read_fasta(IN08 + r'\smbu08\2.Assembly\smbu08.Chromosome.fasta')[0][1]
pl08_recs = read_fasta(IN08 + r'\smbu08\2.Assembly\smbu08.Plasmid.fasta')
pl2 = [r for r in pl08_recs if len(r[1]) == 77437][0][1]
pl1_08 = [r for r in pl08_recs if len(r[1]) == 128837][0][1]
pl1_06 = read_fasta(IN06 + r'\smbu06\2.Assembly\smbu06.Plasmid.fasta')[0][1]
LOG.write('chr06=%d chr08=%d pl2=%d pl1_08=%d pl1_06=%d\n' % (len(chr06), len(chr08), len(pl2), len(pl1_08), len(pl1_06)))
assert len(chr06) - len(chr08) == 38719, 'chr length difference is not 38,719'

# ---------- 1) 染色体比对：从序列重新推导差异区与 att 结构 ----------
write_fasta(os.path.join(REFDIR, '_tmp_chr08.fa'), [('chr08', chr08)])
write_fasta(os.path.join(REFDIR, '_tmp_chr06.fa'), [('chr06', chr06)])
run([BLAST + r'\makeblastdb.exe', '-in', os.path.join(REFDIR, '_tmp_chr06.fa'), '-dbtype', 'nucl',
     '-out', os.path.join(REFDIR, '_tmp_chr06_db')], log=LOG, check=True)
r = run([BLAST + r'\blastn.exe', '-task', 'megablast', '-query', os.path.join(REFDIR, '_tmp_chr08.fa'),
         '-db', os.path.join(REFDIR, '_tmp_chr06_db'),
         '-outfmt', '6 pident length qstart qend sstart send bitscore',
         '-evalue', '1e-5', '-max_hsps', '2000', '-num_threads', '8'], log=LOG, check=True)
big = []
for line in r.stdout.strip().split('\n'):
    if not line:
        continue
    p = line.split('\t')
    if int(p[1]) >= 50000 and int(p[4]) < int(p[5]):
        big.append(dict(pid=float(p[0]), alen=int(p[1]), qs=int(p[2]), qe=int(p[3]),
                        ss=int(p[4]), se=int(p[5]), bits=float(p[6])))
big.sort(key=lambda h: h['qs'])
LOG.write('big HSPs: %s\n' % [(h['qs'], h['qe'], h['ss'], h['se'], h['alen']) for h in big])
assert len(big) == 2, 'expected exactly two backbone HSPs'
front, back = big[0], big[1]
# att 窗口 = 两条 HSP 在 query 坐标上的重叠区（由比对本身导出）
G = back['qs']          # 后段起始（1-based）
F = front['qe']         # 前段终止（1-based）
att_len = F - G + 1
attL = chr06[G - 1:F]
attJ = chr08[G - 1:F]
d = back['ss'] - 1      # 后段在 chr06 的起始（0-based）
attR = chr06[d:d + att_len]
X = chr06[F:d]
element06 = chr06[G - 1:d + att_len]
LOG.write('att_len=%d  attL=%s..%s  X=%d bp  attR=%s..%s  element=%d bp\n' % (
    att_len, G - 1, F, len(X), d, d + att_len, len(element06)))
# 证据表
def diff_positions(a, b):
    n = min(len(a), len(b))
    return [i for i in range(n) if a[i] != b[i]]
dl = diff_positions(attL, attR)
dj_l = diff_positions(attL, attJ)
dj_r = diff_positions(attJ, attR)
with open(os.path.join(ROOT, '01_references', 'att_structure_evidence.tsv'), 'w', encoding='utf-8') as f:
    f.write('item\tvalue\tdetail\n')
    f.write('chr06_len\t%d\t\n' % len(chr06))
    f.write('chr08_len\t%d\t\n' % len(chr08))
    f.write('length_difference\t%d\t与缺失元件长度一致\n' % (len(chr06) - len(chr08)))
    f.write('front_HSP(chr08 vs chr06)\tq%d-%d s%d-%d\t%.4f%% identity, %d bp\n' % (front['qs'], front['qe'], front['ss'], front['se'], front['pid'], front['alen']))
    f.write('back_HSP(chr08 vs chr06)\tq%d-%d s%d-%d\t%.4f%% identity, %d bp\n' % (back['qs'], back['qe'], back['ss'], back['se'], back['pid'], back['alen']))
    f.write('att_window_len(chr08/06 共享)\t%d\t由两条 HSP 的 query 重叠区导出\n' % att_len)
    f.write('attL_chr06\t%d-%d\t%s\n' % (G, F, attL))
    f.write('attJ_chr08\t%d-%d\t%s\n' % (G, F, attJ))
    f.write('attR_chr06\t%d-%d\t%s\n' % (d + 1, d + att_len, attR))
    f.write('attL_vs_attR_diffs\t%d\tpositions(1-based): %s\n' % (len(dl), ';'.join(str(i + 1) for i in dl[:80])))
    f.write('attL_vs_attJ_diffs\t%d\tpositions: %s\n' % (len(dj_l), ';'.join(str(i + 1) for i in dj_l[:80])))
    f.write('attJ_vs_attR_diffs\t%d\tpositions: %s\n' % (len(dj_r), ';'.join(str(i + 1) for i in dj_r[:80])))
    f.write('attJ_byte_equal_attR\t%s\t直接字节比较\n' % (attJ == attR))
    f.write('attJ_byte_equal_attL\t%s\t直接字节比较\n' % (attJ == attL))
    f.write('X_internal_len\t%d\tchr06 %d-%d\n' % (len(X), F + 1, d))
    f.write('element_total_len\t%d\tchr06 %d-%d (attL+X+attR)\n' % (len(element06), G, d + att_len))

# ---------- 2) pl2 结构：单元边界与 junction（自比对自动推导） ----------
write_fasta(os.path.join(REFDIR, '_tmp_pl2.fa'), [('pl2', pl2)])
run([BLAST + r'\makeblastdb.exe', '-in', os.path.join(REFDIR, '_tmp_pl2.fa'), '-dbtype', 'nucl',
     '-out', os.path.join(REFDIR, '_tmp_pl2_db')], log=LOG, check=True)
r = run([BLAST + r'\blastn.exe', '-task', 'dc-megablast', '-query', os.path.join(REFDIR, '_tmp_pl2.fa'),
         '-db', os.path.join(REFDIR, '_tmp_pl2_db'),
         '-outfmt', '6 pident length qstart qend sstart send bitscore qseq sseq',
         '-evalue', '1e-5', '-max_hsps', '20', '-num_threads', '8'], log=LOG, check=True)
h1 = h2 = None   # h1: q 后半 vs s 前半；h2: 对称
for line in r.stdout.strip().split('\n'):
    if not line:
        continue
    p = line.split('\t')
    qs, qe, ss, se = int(p[2]), int(p[3]), int(p[4]), int(p[5])
    if qs > ss and qs > 30000:
        h1 = dict(pid=float(p[0]), alen=int(p[1]), qs=qs, qe=qe, ss=ss, se=se)
    if qs < ss and ss > 30000:
        h2 = dict(pid=float(p[0]), alen=int(p[1]), qs=qs, qe=qe, ss=ss, se=se)
LOG.write('pl2 self HSP: h1=%s h2=%s\n' % (h1, h2))
B = h1['qs'] - 1 if h1 else None      # 单元边界（0-based，unit2 在 pl2 中的起点；h1: q38720-77437 ↔ s1-38719）
unit1 = pl2[:B]
unit2 = pl2[B:]
# 单位 2 相对于单位 1 的差异（逐位；找 1bp 缺失）
mismatch = []
first = None
for i in range(min(len(unit1), len(unit2))):
    if unit1[i] != unit2[i]:
        if first is None:
            first = i
        mismatch.append(i)
        if len(mismatch) > 5:
            break
one_del = False
del_pos = ''
if first is not None:
    k = first
    one_del = (unit2[k:] == unit1[k + 1:len(unit2) + 1]) and unit1[k] == 'A'
    del_pos = k + 1
# 单元与染色体元件的关系（旋转等价）
u_target = revcomp(element06)     # attL+X+attR 的 rc；单元素 U 的候选
rot_eq_u1 = min_rotation(unit1) == min_rotation(revcomp(X + attR))
rot_eq_u1b = min_rotation(unit1) == min_rotation(revcomp(attL + X))
with open(os.path.join(ROOT, '01_references', 'units_and_junctions.tsv'), 'w', encoding='utf-8') as f:
    f.write('item\tvalue\tdetail\n')
    f.write('pl2_len\t%d\t\n' % len(pl2))
    f.write('self_HSP_first_half_unit1_to_unit2\tq%d-%d s%d-%d %.4f%%\t融合单元取向\n' % (h1['qs'], h1['qe'], h1['ss'], h1['se'], h1['pid']))
    f.write('self_HSP_second\tq%d-%d s%d-%d %.4f%%\t\n' % (h2['qs'], h2['qe'], h2['ss'], h2['se'], h2['pid']))
    f.write('unit_boundary_0based\t%d\tunit1=[0,%d), unit2=[%d,%d)\n' % (B, B, B, len(pl2)))
    f.write('unit1_len\t%d\t\n' % len(unit1))
    f.write('unit2_len\t%d\t\n' % len(unit2))
    f.write('unit2_minus_unit1\t1bp deletion verified=%s at unit1 position %s (base=%s)\n' % (one_del, del_pos, unit1[del_pos - 1] if del_pos else ''))
    f.write('unit1_rotation_eq_revcomp(attL+X)\t%s\t（U = attL+X 的反向互补）\n' % rot_eq_u1b)
    f.write('unit1_rotation_eq_revcomp(X+attR)\t%s\t\n' % rot_eq_u1)
    f.write('unit2_to_unit1_boundary(origin)\t0/pl2_len\t环状回接连接位点\n')

# ---------- 3) 参考面板构建 ----------
FL = 5000
j_att = (G - 1 + F) // 2                      # attJ 中点（0-based, chr08）
attJ_junction = chr08[j_att - FL:j_att + FL]
j_fwd = B                                     # unit1→unit2 边界（pl2 0-based）
unit1_to_unit2 = pl2[j_fwd - FL:j_fwd + FL]
back_junction = pl2[-FL:] + pl2[:FL]          # unit2→unit1（环状原点）连接位点
retained_window = chr06[max(0, G - 1 - FL):d + att_len + FL]
# 单体与二聚体模型（环状以两处线性化点各有一条；线性=组装原样）
monomer_A = unit1
monomer_B = unit1[len(unit1) // 2:] + unit1[:len(unit1) // 2]
dimer_asassembled = pl2                        # unit1→unit2 junction 内部（原点在回接处）
dimer_back_internal = pl2[B:] + pl2[:B]        # unit2→unit1 junction 内部
# 阴性对照 junction（GC 匹配、随机位置，避开已知 junction ±100 kb）
import random
rng = random.Random(20261002)
gc_target = (attJ_junction.count('G') + attJ_junction.count('C')) / len(attJ_junction)
neg = []
tries = 0
while len(neg) < 3 and tries < 20000:
    tries += 1
    pos = rng.randrange(FL, len(chr08) - FL)
    if abs(pos - j_att) < 100000:
        continue
    w = chr08[pos - FL:pos + FL]
    gc = (w.count('G') + w.count('C')) / len(w)
    if abs(gc - gc_target) <= 0.005:
        neg.append(('neg_ctrl_chr08_%d' % (len(neg) + 1), w, pos))
# 单核苷酸打乱的 attJ 窗口（保守阴性对照）
ww = list(attJ_junction)
rng.shuffle(ww)
neg_shuf = ('neg_ctrl_shuffled_attJ', ''.join(ww), -1)

refs = [
    ('chr08_full', chr08, 'smbu08 chromosome 全序列（含 attJ 状态）'),
    ('chr06_full', chr06, 'smbu06 chromosome 全序列（保留态 attL-X-attR）'),
    ('pl2_asassembled', pl2, 'smbu08 77,437 bp 复制子，组装原样（unit1→unit2 junction 在 %d 内部；回接连接在 0/end）' % j_fwd),
    ('pl1_smbu08', pl1_08, '128,837 bp 共享质粒（无关对照）'),
    ('attJ_junction', attJ_junction, 'chr08 以 attJ 中点为中心 ±%d bp；junction 位于 %d' % (FL, FL)),
    ('retained_window', retained_window, 'chr06 保留态窗口：A尾%dbp + attL+X+attR + B头%dbp；junction 覆盖 %d-%d' % (FL, FL, FL, FL + att_len)),
    ('unit1', unit1, 'pl2 单元 1（0–%d）' % B),
    ('unit2', unit2, 'pl2 单元 2（%d–%d）' % (B, len(pl2))),
    ('unit1_to_unit2_junction', unit1_to_unit2, 'pl2 上 unit1→unit2 边界 ±%d bp；junction 位于 %d' % (FL, FL)),
    ('unit2_to_unit1_circular_junction', back_junction, 'pl2 环状原点（unit2 尾 + unit1 头）±%d bp；junction 位于 %d' % (FL, FL)),
    ('monomer_circle_repA', monomer_A, '单体环状模型（线性化点=原点）'),
    ('monomer_circle_repB', monomer_B, '单体环状模型（线性化点=半程旋转）'),
    ('dimer_back_internal', dimer_back_internal, '二聚体环状模型（线性化于 unit1→unit2 边界，使回接连接点内部化）'),
    ('neg_attJ_shuffled', neg_shuf[1], 'attJ 窗口单核苷酸打乱（阴性对照；junction=%d）' % FL),
]
for name, seq, note in neg:
    refs.append((name, seq, 'GC 匹配随机阴性对照（chr08 pos %d，GC=%.4f；junction=%d）' % (note, (seq.count('G') + seq.count('C')) / len(seq), FL)))

manifest = []
for name, seq, note in refs:
    p = os.path.join(REFDIR, name + '.fasta')
    write_fasta(p, [(name, seq)])
    h = hashlib.sha256((seq + '\n').encode()).hexdigest()
    manifest.append((name, len(seq), hashlib.sha256(open(p,'rb').read()).hexdigest(), note))
    LOG.write('ref %s len=%d\n' % (name, len(seq)))

with open(os.path.join(ROOT, '01_references', 'reference_construction.tsv'), 'w', encoding='utf-8') as f:
    f.write('reference\tlength_bp\tfasta_sha256\tconstruction\n')
    for r_ in manifest:
        f.write('\t'.join(str(x) for x in r_) + '\n')

# junction 定义表（供映射分析使用）
jrows = [
    ('attJ_junction', 'attJ_chr08_state', FL, 'chr08 attJ 中点；两侧 %dbp 侧翼' % FL),
    ('unit1_to_unit2_junction', 'dimer_forward', FL, 'pl2 unit1→unit2 边界；两侧 %dbp' % FL),
    ('unit2_to_unit1_circular_junction', 'dimer_back_wrap', FL, 'pl2 原点回接；两侧 %dbp' % FL),
    ('retained_window', 'chr06_retained_state', FL + att_len // 2, 'chr06 保留元件窗口（attL-X-attR）'),
    ('neg_attJ_shuffled', 'negative_control', FL, '打乱序列阴性对照'),
]
for name, seq, note in neg:
    jrows.append((name, 'negative_control', FL, 'GC 匹配随机对照'))
with open(os.path.join(ROOT, '06_tables', 'junction_definitions.tsv'), 'w', encoding='utf-8') as f:
    f.write('reference\tjunction_type\tjunction_pos(0-based)\tnote\n')
    for r_ in jrows:
        f.write('\t'.join(str(x) for x in r_) + '\n')

for tmp in ['_tmp_chr08.fa', '_tmp_chr06.fa', '_tmp_pl2.fa']:
    for ext in ['', '.nhr', '.nin', '.nsq', '.ndb', '.not', '.ntf', '.nto', '.njs', '.nhd', '.nhi', '.nnd', '.nsi', '.nog', '.nsd']:
        p = os.path.join(REFDIR, tmp + ext if ext else tmp)
        if os.path.exists(p):
            os.remove(p)
LOG.write('DONE B=%d att_len=%d j_att=%d\n' % (B, att_len, j_att))
print('build_references done: B=%d att_len=%d unit2_1bpdel=%s attJ==attR:%s' % (B, att_len, one_del, attJ == attR))
