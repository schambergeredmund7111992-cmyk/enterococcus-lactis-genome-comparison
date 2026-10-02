# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
chr_junction.py — 染色体差异的字节级精确分解（模块 03 核心）
结构（全部字节级验证）：
  smbu06 chr = A + attL(146) + X'(38,573) + attR(146) + B
  smbu08 chr = A + attJ(146) + B,  attJ == attR（字节相等）
  → smbu08 缺失 attL+X'（38,719 bp），另有 1 个远端 SNP
输出：
  03_replicon_compare/chr_difference_summary.tsv
  03_replicon_compare/chr_snp_list.tsv
  03_replicon_compare/att_region_analysis.tsv
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
import numpy as np
from p3lib import WORK, ASM06, ASM08, read_fasta, revcomp, open_log

OUT = os.path.join(WORK, '03_replicon_compare')
LOG = open_log(os.path.join(OUT, 'logs', 'chr_junction.log'))

s = read_fasta(ASM06['Chromosome1'])[0][1]
q = read_fasta(ASM08['Chromosome1'])[0][1]
LOG.write('s len %d, q len %d, diff %d\n' % (len(s), len(q), len(s) - len(q)))

# ---------- 1) 精确分解 ----------
# 从比对已知：q[0:1255782]==s[0:1255782]; q[1255782:1255928]==s[1294501:1294647];
# q[1255928:]==s[1294647:]. 逐段字节验证。
Aq, As = q[0:1255782], s[0:1255782]
attJ, attR = q[1255782:1255928], s[1294501:1294647]
Bq, Bs = q[1255928:], s[1294647:]
checks = [
    ('A 区段相等 (q[1..1255782]==s[1..1255782])', Aq == As),
    ('attJ == attR (q[1255783..1255928]==s[1294502..1294647])', attJ == attR),
    ('B 区段相等 (q[1255929..]==s[1294648..])', Bq == Bs),
]
attL = s[1255782:1255928]
Xp = s[1255928:1294501]
LOG.write('attL len %d, Xp len %d, attR len %d\n' % (len(attL), len(Xp), len(attR)))
checks.append(('长度守恒: len(attL)+len(Xp)+len(attR) == 38719', len(attL) + len(Xp) + len(attR) == 38719))
for name, ok in checks:
    LOG.write('CHECK %s : %s\n' % (name, ok))

# ---------- 2) B 区段内的 SNP（逐字节） ----------
aq = np.frombuffer(Bq.encode(), dtype=np.uint8)
bs = np.frombuffer(Bs.encode(), dtype=np.uint8)
assert len(aq) == len(bs)
mm = np.nonzero(aq != bs)[0]
snp_rows = []
for i in mm:
    qpos = 1255928 + int(i) + 1     # 1-based
    spos = 1294647 + int(i) + 1
    snp_rows.append((qpos, spos, chr(aq[i]), chr(bs[i])))
LOG.write('B 区段 SNP 数: %d -> %s\n' % (len(mm), snp_rows))

# A 区段全面扫描（应为 0）
aa = np.frombuffer(Aq.encode(), dtype=np.uint8)
bb = np.frombuffer(As.encode(), dtype=np.uint8)
mmA = np.nonzero(aa != bb)[0]
LOG.write('A 区段 SNP 数: %d\n' % len(mmA))

# ---------- 3) att 区分析 ----------
# attL vs attR 差异
diff_att = [(i + 1, attL[i], attR[i]) for i in range(146) if attL[i] != attR[i]]
# 最长公共子串（attL vs attR）
def lcs(a, b):
    best = ''
    for L in range(min(len(a), len(b)), 9, -1):
        found = False
        for i in range(len(a) - L + 1):
            if a[i:i + L] in b:
                best = a[i:i + L]
                found = True
                break
        if found:
            return L, best
    return 0, ''
Lc, core = lcs(attL, attR)
LOG.write('attL vs attR: 差异碱基 %d/146; 最长公共子串 %dbp: %s\n' % (len(diff_att), Lc, core))

with open(os.path.join(OUT, 'chr_difference_summary.tsv'), 'w', encoding='utf-8') as f:
    f.write('item\tvalue\tdetail\n')
    f.write('chromosome_len_smbu06\t%d\t\n' % len(s))
    f.write('chromosome_len_smbu08\t%d\t\n' % len(q))
    f.write('length_difference_bp\t%d\t单缺失，与attL+X\'长度一致\n' % (len(s) - len(q)))
    f.write('prophage_segment_in_smbu06\t1,255,783-1,294,501\t38,719 bp = attL(146)+X\'(38,573)，smbu08缺失\n')
    f.write('prophage_segment_alt_frame\t1,255,929-1,294,647\t38,719 bp = X\'+attR（等价框，差一个att拷贝的归属）\n')
    f.write('left_att_in_smbu06(attL)\ts[1255782:1255928]\t146 bp\n')
    f.write('right_att_in_smbu06(attR)\ts[1294501:1294647]\t146 bp，与smbu08连接序列字节相等\n')
    f.write('smbu08_junction\tq[1255782:1255928]\t146 bp，== attR（字节相等）\n')
    f.write('attL_vs_attR_diffs\t%d/146\t见 att_region_analysis.tsv\n' % len(diff_att))
    f.write('att_core_lcs\t%d bp\t%s\n' % (Lc, core))
    f.write('snps_outside_prophage\t%d\t见 chr_snp_list.tsv\n' % len(mm))
    f.write('A_region_snps\t%d\t\n' % len(mmA))
    f.write('B_region_snps\t%d\t\n' % len(mm))

with open(os.path.join(OUT, 'chr_snp_list.tsv'), 'w', encoding='utf-8') as f:
    f.write('smbu08_chr_pos(1-based)\tsmbu06_chr_pos(1-based)\tsmbu08_base\tsmbu06_base\n')
    for r in snp_rows:
        f.write('%d\t%d\t%s\t%s\n' % r)

with open(os.path.join(OUT, 'att_region_analysis.tsv'), 'w', encoding='utf-8') as f:
    f.write('pos_in_att(1-based)\tattL_base\tattR_base\n')
    for r in diff_att:
        f.write('%d\t%s\t%s\n' % r)
    # 附序列全文，便于复核
    f.write('#attL_seq\t%s\n' % attL)
    f.write('#attR_seq\t%s\n' % attR)
    f.write('#attJ_seq(==attR)\t%s\n' % attJ)

LOG.write('DONE\n')
print('chr_junction done. SNP counts: A=%d B=%d att_diffs=%d' % (len(mmA), len(mm), len(diff_att)))
