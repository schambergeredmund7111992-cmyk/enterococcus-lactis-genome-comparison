# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
pl2_origin.py — smbu08 Plasmid2 (77,437 bp) 结构鉴定（模块 03/04 交汇）
验证：
  1) plasmid2 = U 的串联二聚体（U = revcomp(attL+X')，38,719 bp，两者旋转等价）
  2) 第二拷贝相对第一拷贝缺失 1 bp（对齐级证据）
  3) plasmid2 对 smbu06 染色体的匹配区间与 att 结构
输出：03_replicon_compare/pl2_structure_summary.tsv
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
import numpy as np
from p3lib import (WORK, ASM06, ASM08, BLASTN, read_fasta, revcomp,
                   min_rotation, open_log)

OUT = os.path.join(WORK, '03_replicon_compare')
LOG = open_log(os.path.join(OUT, 'logs', 'pl2_origin.log'))

s = read_fasta(ASM06['Chromosome1'])[0][1]
pl2 = [r for r in read_fasta(ASM08['Plasmid1']) if 'Plasmid2' in r[0]][0][1]
h1, h2 = pl2[:38719], pl2[38719:]
U = revcomp(s[1255782:1294501])          # attL + X' 的反向互补 = 38,719

rows = []
rot_h1_eq = min_rotation(h1) == min_rotation(U)
rot_h2_eq = min_rotation(h2) == min_rotation(U)
# 找 h1 相对 U 的旋转量（O(n)，用双串 find）
rot = (U + U).find(h1) if len(h1) == len(U) else None
rows.append(('plasmid2_length', len(pl2)))
rows.append(('unit_length', len(U)))
rows.append(('half1_equals_rotation_of_U', rot_h1_eq))
rows.append(('half1_rotation_offset', rot))
rows.append(('half2_equals_rotation_of_U_minus_1bp', rot_h2_eq))
rows.append(('dimer_junction_1bp_deletion', not rot_h2_eq and len(h2) == len(h1) - 1))
# 二聚体第二拷贝缺 1 bp 的直接证据（numpy，O(n)）：
# 找 h2 与 h1 的第一个错位点 p，验证 h2[p:] == h1[p+1:len(h2)+1]
n = len(h2)
aa = np.frombuffer(h2.encode(), dtype=np.uint8)
bb = np.frombuffer(h1[:n].encode(), dtype=np.uint8)
p = int(np.nonzero(aa != bb)[0][0])
one_del = (h2[p:] == h1[p + 1:n + 1]) and (h1[p] == 'A')
rows.append(('single_deletion_verified', one_del))
rows.append(('deletion_position_in_unit(1-based)', p + 1))
rows.append(('deleted_base', h1[p]))

# plasmid2 与 smbu06 染色体区域的覆盖
import subprocess
pl2fa = os.path.join(OUT, 'alignments', 'smbu08_pl2.fa')
db = os.path.join(OUT, 'alignments', 'smbu06_chr_db')
r = subprocess.run([BLASTN, '-task', 'dc-megablast', '-query', pl2fa, '-db', db,
                    '-outfmt', '6 qstart qend sstart send length pident mismatch gaps',
                    '-evalue', '1e-5', '-max_hsps', '50'], capture_output=True, text=True)
big = []
for line in r.stdout.strip().split('\n'):
    if not line:
        continue
    p = line.split('\t')
    if int(p[4]) >= 1000:
        big.append(p)
LOG.write('big HSPs (>=1kb):\n')
for p in big:
    LOG.write('  q%s-%s s%s-%s len=%s pid=%s\n' % tuple(p[:6]))
rows.append(('pl2_vs_smbu06chr_bigHSPs', ';'.join('q%s-%s:s%s-%s(%s%%,%sbp)' % (p[0], p[1], p[2], p[3], p[5], p[4]) for p in big)))

with open(os.path.join(OUT, 'pl2_structure_summary.tsv'), 'w', encoding='utf-8') as f:
    f.write('item\tvalue\n')
    for k, v in rows:
        f.write('%s\t%s\n' % (k, v))

LOG.write('\n'.join('%s = %s' % (k, v) for k, v in rows) + '\nDONE\n')
print('pl2_origin done:', rows[:6])
