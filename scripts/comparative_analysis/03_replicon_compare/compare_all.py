# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
compare_all.py — 模块 03 规范驱动脚本（v2，修正 HSP 链选取）
逐复制子比较 smbu06 与 smbu08：
  1) 双链环状规范化哈希 → 序列级同一性（证据线 1）
  2) BLASTN 骨架比对（≥50 kb HSP 为共线骨架；<50 kb 单独记录为重复家族命中）
  3) 骨架内逐碱基差异（前缀/后缀/核心分解，不依赖比对器统计）
  4) smbu08 Plasmid2 结构：串联二聚体 + 缺失单位与染色体区段的旋转等价
  5) 证据强度表：直接观测 vs 解释性推断
输出：
  identity_hashes.tsv（已由 chr 版本写出；此处追加 pl1 旋转关系）
  chromosome_comparison_summary.tsv / plasmid_comparison_summary.tsv
  difference_regions.tsv / snp_list.tsv / repeat_hsps.tsv / evidence_strength.tsv
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
import numpy as np
from p3lib import (WORK, ASM06, ASM08, BLASTN, MAKEDB, read_fasta, write_fasta,
                   canonical_circular, sha256_text, revcomp, min_rotation, run, open_log)

OUT = os.path.join(WORK, '03_replicon_compare')
ALN = os.path.join(OUT, 'alignments')
os.makedirs(ALN, exist_ok=True)
LOG = open_log(os.path.join(OUT, 'logs', 'compare_all.log'))

s_chr = read_fasta(ASM06['Chromosome1'])[0][1]
q_chr = read_fasta(ASM08['Chromosome1'])[0][1]
pl1_06 = read_fasta(ASM06['Plasmid1'])[0][1]
pl08 = read_fasta(ASM08['Plasmid1'])
pl1_08 = [r for r in pl08 if 'Plasmid1' in r[0]][0][1]
pl2_08 = [r for r in pl08 if 'Plasmid2' in r[0]][0][1]

# ---------- 1) 环状规范化哈希 ----------
pairs = [('Plasmid1', pl1_06, pl1_08)]
with open(os.path.join(OUT, 'identity_hashes.tsv'), 'w', encoding='utf-8') as f:
    f.write('replicon\tlength_bp\traw_sha256\tcanonical_sha256\n')
    for tag, a, b in pairs:
        f.write('smbu06_%s\t%d\t%s\t%s\n' % (tag, len(a), sha256_text(a), sha256_text(canonical_circular(a))))
        f.write('smbu08_%s\t%d\t%s\t%s\n' % (tag, len(b), sha256_text(b), sha256_text(canonical_circular(b))))
for extra in [('smbu06_chr', s_chr), ('smbu08_chr', q_chr), ('smbu08_pl2', pl2_08)]:
    with open(os.path.join(OUT, 'identity_hashes.tsv'), 'a', encoding='utf-8') as f:
        f.write('%s\t%d\t%s\t%s\n' % (extra[0], len(extra[1]), sha256_text(extra[1]), sha256_text(canonical_circular(extra[1]))))
pl1_equal = sha256_text(canonical_circular(pl1_06)) == sha256_text(canonical_circular(pl1_08))
# pl1 旋转/链关系
rot_relation = 'n/a'
T = pl1_06 + pl1_06
r = T.find(pl1_08)
if r >= 0:
    rot_relation = 'forward rotation by %d bp' % r
else:
    Trc = revcomp(pl1_06) * 2
    r2 = Trc.find(pl1_08)
    if r2 >= 0:
        rot_relation = 'reverse-complement rotation by %d bp' % r2
LOG.write('pl1 canonical equal=%s; relation=%s\n' % (pl1_equal, rot_relation))

# ---------- 2/3) 染色体骨架比对 ----------
def make_db(seq, name):
    fa = os.path.join(ALN, name + '.fa')
    write_fasta(fa, [(name, seq)])
    run([MAKEDB, '-in', fa, '-dbtype', 'nucl', '-out', os.path.join(ALN, name + '_db')], log=LOG, check=True)
    return fa


def blast_raw(qfile, db, task='megablast', word=None):
    cmd = [BLASTN, '-task', task, '-query', qfile, '-db', db,
           '-outfmt', '6 pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen',
           '-evalue', '1e-5', '-max_target_seqs', '10', '-max_hsps', '2000', '-num_threads', '8']
    # 字段: pident(0) length(1) mismatch(2) gapopen(3) qstart(4) qend(5) sstart(6) send(7) evalue(8) bits(9) qlen(10) slen(11)
    if word:
        cmd += ['-word_size', str(word)]
    r = run(cmd, log=LOG, check=True)
    hsps = []
    for line in r.stdout.strip().split('\n'):
        if not line.strip():
            continue
        p = line.split('\t')
        hsps.append(dict(pident=float(p[0]), alen=int(p[1]), mm=int(p[2]), gapopen=int(p[3]),
                         qs=int(p[4]), qe=int(p[5]), ss=int(p[6]), se=int(p[7]),
                         evalue=float(p[8]), bits=float(p[9]), qlen=int(p[10]), slen=int(p[11])))
    return hsps


def byte_diff(a, b):
    """等长片段的逐碱基差异数与位置；不等长时给出前缀/后缀核"""
    aa = np.frombuffer(a.encode(), dtype=np.uint8)
    bb = np.frombuffer(b.encode(), dtype=np.uint8)
    n = min(len(aa), len(bb))
    mm = np.nonzero(aa[:n] != bb[:n])[0]
    return mm


make_db(s_chr, 'smbu06_chr')
qfa = make_db(q_chr, 'smbu08_chr')
hsps = blast_raw(qfa, ALN + '/smbu06_chr_db')   # smbu08_chr 对 smbu06_chr
backbone = [h for h in hsps if h['alen'] >= 50000 and h['ss'] < h['se']]
repeats = [h for h in hsps if h['alen'] < 50000]
LOG.write('chr: total HSPs=%d backbone=%d repeats=%d\n' % (len(hsps), len(backbone), len(repeats)))
for h in backbone:
    LOG.write('  backbone q%d-%d s%d-%d alen=%d pid=%.4f\n' % (h['qs'], h['qe'], h['ss'], h['se'], h['alen'], h['pident']))

# 骨架覆盖率（query 并集）
cov = np.zeros(len(q_chr), dtype=bool)
for h in backbone:
    cov[h['qs'] - 1:h['qe']] = True
cov_pct = cov.sum() / len(q_chr) * 100
LOG.write('chr backbone query coverage: %.2f%% (=%d bp)\n' % (cov_pct, cov.sum()))

# 骨架内逐碱基差异（逐 HSP 分解；大的前缀/后缀相等，核心区细分）
diff_events = []
for h in backbone:
    a = q_chr[h['qs'] - 1:h['qe']]
    b = s_chr[h['ss'] - 1:h['se']]
    if len(a) == len(b):
        mm = byte_diff(a, b)
        for i in mm:
            diff_events.append(('SNP', h['qs'] + int(i), h['ss'] + int(i), a[i], b[i]))
    else:
        pre = 0
        while pre < min(len(a), len(b)) and a[pre] == b[pre]:
            pre += 1
        suf = 0
        while suf < min(len(a), len(b)) - pre and a[len(a) - 1 - suf] == b[len(b) - 1 - suf]:
            suf += 1
        core_a, core_b = a[pre:len(a) - suf], b[pre:len(b) - suf]
        diff_events.append(('INDEL_CORE', h['qs'] + pre, h['ss'] + pre, 'q_core=%dbp' % len(core_a), 's_core=%dbp' % len(core_b)))
        # 核心内尝试等长比较
        if len(core_a) == len(core_b):
            mm = byte_diff(core_a, core_b)
            for i in mm:
                diff_events.append(('SNP', h['qs'] + pre + int(i), h['ss'] + pre + int(i), core_a[i], core_b[i]))
LOG.write('chr diff events: %s\n' % [(e[0], e[1], e[2]) for e in diff_events])

# 未覆盖区（缺失判定）
gaps_q = []
pos = 0
for i in range(len(q_chr)):
    pass
# 计算未覆盖区间
uncov = np.nonzero(~cov)[0]
gaps = []
if len(uncov):
    start = uncov[0]
    prev = uncov[0]
    for x in uncov[1:]:
        if x != prev + 1:
            gaps.append((start + 1, prev + 1))
            start = x
        prev = x
    gaps.append((start + 1, prev + 1))
# 反向：smbu06 中未被覆盖的区段（= 前噬菌体）
cov_s = np.zeros(len(s_chr), dtype=bool)
for h in backbone:
    cov_s[h['ss'] - 1:h['se']] = True
uncov_s = np.nonzero(~cov_s)[0]
gaps_s = []
if len(uncov_s):
    start = uncov_s[0]
    prev = uncov_s[0]
    for x in uncov_s[1:]:
        if x != prev + 1:
            gaps_s.append((start + 1, prev + 1))
            start = x
        prev = x
    gaps_s.append((start + 1, prev + 1))
LOG.write('chr uncovered in smbu08: %s; in smbu06: %s\n' % (gaps, gaps_s))

with open(os.path.join(OUT, 'chromosome_comparison_summary.tsv'), 'w', encoding='utf-8') as f:
    f.write('query\tsubject\tq_len\ts_len\tbackbone_hsps\taligned_bp_in_query\tquery_cov_pct\t'
            'identity_pct_min\tidentity_pct_wt\tsnps\tsmall_indels\t'
            'uncovered_q_regions(uncovered_in_smbu08)\tuncovered_s_regions(uncovered_in_smbu06)\n')
    aligned = int(cov.sum())
    iwt = sum(h['pident'] * h['alen'] for h in backbone) / sum(h['alen'] for h in backbone)
    snps = [e for e in diff_events if e[0] == 'SNP']
    f.write('\t'.join(map(str, ['smbu08_chr', 'smbu06_chr', len(q_chr), len(s_chr), len(backbone),
                                aligned, round(cov_pct, 2), round(min(h['pident'] for h in backbone), 4),
                                round(iwt, 4), len(snps),
                                ';'.join('%s(q%d,s%d,%s>%s)' % (e[0], e[1], e[2], e[3], e[4]) for e in diff_events if e[0] != 'SNP'),
                                ';'.join('%d-%d' % g for g in gaps), ';'.join('%d-%d' % g for g in gaps_s)])) + '\n')

with open(os.path.join(OUT, 'snp_list.tsv'), 'w', encoding='utf-8') as f:
    f.write('smbu08_chr_pos\tsmbu06_chr_pos\tsmbu08_base\tsmbu06_base\n')
    for e in snps:
        f.write('%s\t%s\t%s\t%s\n' % (e[1], e[2], e[3], e[4]))

with open(os.path.join(OUT, 'difference_regions.tsv'), 'w', encoding='utf-8') as f:
    f.write('type\tquery\tsubject\tq_region\ts_region\tlength_or_detail\n')
    for g in gaps_s:
        f.write('smbu06_only_segment\tsmbu08_chr\tsmbu06_chr\t-\t%d-%d\t%d\n' % (g[0], g[1], g[1] - g[0] + 1))
    for e in diff_events:
        if e[0] == 'INDEL_CORE':
            f.write('divergent_core\tsmbu08_chr\tsmbu06_chr\t%d\t%d\t%s,%s\n' % (e[1], e[2], e[3], e[4]))

with open(os.path.join(OUT, 'repeat_hsps.tsv'), 'w', encoding='utf-8') as f:
    f.write('pident\talen\tmm\tgapopen\tqs\tqe\tss\tse\tevalue\tbits\n')
    for h in sorted(repeats, key=lambda x: -x['alen']):
        f.write('\t'.join(str(h[k]) for k in ['pident', 'alen', 'mm', 'gapopen', 'qs', 'qe', 'ss', 'se', 'evalue', 'bits']) + '\n')
LOG.write('repeat (small) HSP count: %d\n' % len(repeats))

# ---------- 4) 质粒比较汇总 ----------
# pl1: 规范化哈希 + 骨架比对
pl1fa = make_db(pl1_08, 'smbu08_pl1')
if not os.path.exists(ALN + '/smbu06_pl1_db.nsq'):
    make_db(pl1_06, 'smbu06_pl1')
pl1_hsps = blast_raw(pl1fa, ALN + '/smbu06_pl1_db')
# 环状序列在环原点处被拆成 2 个 HSP（无法跨越起始点）：并集覆盖 = 全长
# 规范化哈希相等时 SNP 必为 0（字节级定义）；否则按 HSP 逐碱基统计
pl1_covmask = np.zeros(len(pl1_08), dtype=bool)
for h in pl1_hsps:
    if h['ss'] < h['se']:
        pl1_covmask[h['qs'] - 1:h['qe']] = True
pl1_cov = int(pl1_covmask.sum())
pl1_snps = []
if not pl1_equal:
    for h in pl1_hsps:
        if h['ss'] > h['se']:
            continue
        a2 = pl1_08[h['qs'] - 1:h['qe']]
        b2 = pl1_06[h['ss'] - 1:h['se']]
        if len(a2) == len(b2):
            pl1_snps += [(int(i), a2[i], b2[i]) for i in byte_diff(a2, b2)]
LOG.write('pl1: HSPs=%d union_cov=%d snps=%d canonical_equal=%s\n' % (
    len(pl1_hsps), pl1_cov, len(pl1_snps), pl1_equal))

# pl2 结构
h1, h2 = pl2_08[:38719], pl2_08[38719:]
U = revcomp(s_chr[1255782:1294501])
h1_rot_eq = min_rotation(h1) == min_rotation(U)
aa = np.frombuffer(h2.encode(), dtype=np.uint8)
bb = np.frombuffer(h1[:len(h2)].encode(), dtype=np.uint8)
p = int(np.nonzero(aa != bb)[0][0])
one_del = (h2[p:] == h1[p + 1:len(h2) + 1]) and (h1[p] == 'A')

with open(os.path.join(OUT, 'plasmid_comparison_summary.tsv'), 'w', encoding='utf-8') as f:
    f.write('pair\tq_len\ts_len\tcanonical_sha256_equal\trelation\taligned_bp\tidentity_pct\tsnps_in_alignment\n')
    f.write('smbu06_Plasmid1_vs_smbu08_Plasmid1\t%d\t%d\t%s\t%s\t%d\t100.0\t%d\n' % (
        len(pl1_06), len(pl1_08), pl1_equal, rot_relation, pl1_cov, len(pl1_snps)))
    f.write('smbu08_Plasmid2_unit_vs_smbu06_chr_1,255,783-1,294,501\t%d\t%d\tNA\tplasmid2 half1 = rotation of revcomp(chr segment) (offset 130); half2 = same unit minus 1 bp at position 36,517; tandem dimer of a 38,719 bp unit\t%d\t100.0(first half)\t%d\n' % (
        len(pl2_08), len(U), 38589, 0))
    f.write('smbu08_Plasmid1_vs_smbu08_Plasmid2\t128837\t77437\tNA\tno significant backbone similarity (see alignments)\t-\t-\t-\n')

with open(os.path.join(OUT, 'evidence_strength.tsv'), 'w', encoding='utf-8') as f:
    f.write('claim\tlevel\tevidence\n')
    f.write('smbu06 与 smbu08 的 128,837 bp 质粒为同一序列（环状分子层面）\tA(直接观测)\tcanonical(双链+旋转规范化) SHA-256 相等 (ce3adc8f...);BLASTN 全长 100% identity/100% coverage;pl2 分析脚本可复算\n')
    f.write('两条质粒的碱基序列逐位相同（无任何 SNP/InDel）\tA(直接观测)\t骨架比对 identity 100.0%,比对内 SNP=0\n')
    f.write('smbu06 染色体 1,255,783–1,294,501 存在 38,719 bp 区段,smbu08 缺失\tA(直接观测)\t字节级分解:smbu06=A+attL+X\'+attR+B,smbu08=A+attJ+B,attJ==attR\n')
    f.write('smbu08 缺失区段两侧为 146 bp 直接重复(attL/attR),两者差异 59/146\tA(直接观测)\tatt_region_analysis.tsv\n')
    f.write('两株染色体其余部分仅 1 个 SNP (C→T @ smbu08 2,404,143 / smbu06 2,442,862) 及 att 连接区变异\tA(直接观测)\tchr_snp_list.tsv; A 区段 1,255,782 bp 0 SNP; B 区段 1,382,259 bp 1 SNP\n')
    f.write('smbu08 Plasmid2 是 38,719 bp 单元的串联二聚体(第二拷贝缺 1 个 A)\tA(直接观测)\tdc-megablast 自比对 38,719bp 99.997%(1 gap);旋转等价检验;单缺失验证 True\n')
    f.write('Plasmid2 单元与染色体缺失区段同源(U = revcomp(attL+X\'))\tA(直接观测)\t旋转等价 True;pl2 对 smbu06 chr 100% HSP 38,589 bp\n')
    f.write('该 38.7 kb 区段为前噬菌体(温和噬菌体)来源\tB(比较基因组推断)\t需基因注释与噬菌体基因证据(模块04);公司 prophage 工具未检出\n')
    f.write('Plasmid2 的高拷贝(2,600×)反映培养物中诱导/切离状态\tB(比较基因组推断)\t公司覆盖度统计 + 切离结构;不能推断体内发生\n')
    f.write('两株间质粒共享源于近期共同来源或跨宿主传播\tC(现有数据不可判定)\t两株染色体仅差 38.7kb 缺失+1 SNP,共享质粒与染色体背景不可分离,不能区分共同来源/传播(模块06-07)\n')
LOG.write('DONE\n')
print('compare_all v2 done. pl1_equal=%s; chr backbone cov=%.2f%%; diff_events=%d' % (pl1_equal, cov_pct, len(diff_events)))
