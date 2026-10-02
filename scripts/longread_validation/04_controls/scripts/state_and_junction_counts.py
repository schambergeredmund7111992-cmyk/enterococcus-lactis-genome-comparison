# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
state_and_junction_counts.py — 染色体状态与 junction 支持的统一计数（PAF 路径）
对 smbu08 与 smbu06（阳性对照）分别统计：
  保留态接合  chr06 flank|attL @ S_pos（仅保留态模板含此邻接）
  切离态接合  chr08 flank|attJ @ S_pos（仅切离态模板含此邻接）
  保留态右接合 chr06 X|attR @ S_pos+B（仅保留态模板）
  切离态右接合 chr08 attJ|flank @ S_pos+A（仅切离态模板）
  质粒正向 junction  pl2 u1|u2 @ 38719
  质粒回接 junction  dimer_back_internal @ 38718（旋转帧）
判据：单条比对连续跨越（M 覆盖 junction），两侧匹配锚定 ≥1kb（严格 ≥3kb），mapq≥20。
输出：04_controls/state_junction_counts.tsv
"""
import gzip
import os
import re
from collections import defaultdict, Counter

ATTJ = r'<PROJECT_ROOT>\longread_validation'
CIG_RE = re.compile(r'(\d+)([MIDNSHP=X])')
S_pos = 1255782
B = 38719
A = 146

JOINTS = {
    'retained_left(chr06 flank|attL)': ('chr06_full', S_pos),
    'excised_left(chr08 flank|attJ)': ('chr08_full', S_pos),
    'retained_right(chr06 X|attR)': ('chr06_full', S_pos + B),
    'excised_right(chr08 attJ|flank)': ('chr08_full', S_pos + A),
    'plasmid_fwd(u1|u2)': ('pl2_asassembled', 38719),
    'plasmid_back(origin)': ('dimer_back_internal', 38718),
}


def anchors(ops, r0, jpos):
    pos = r0
    l = r = 0
    found = False
    for n, o in ops:
        if o in 'M=X':
            if pos <= jpos < pos + n:
                found = True
                l += jpos - pos
                r += pos + n - jpos - 1
            elif pos + n <= jpos:
                l += n
            else:
                r += n
            pos += n
        elif o in 'DN':
            if pos <= jpos < pos + n:
                return None
            pos += n
    return (l, r) if found else None


DATASETS = [
    ('smbu08', os.path.join(ATTJ, '02_mapping', 'full_db.paf.gz')),
    ('smbu06_control', os.path.join(ATTJ, '04_controls', 'smbu06_db.paf.gz')),
    ('smbu08_pl2frames', os.path.join(ATTJ, '04_controls', 'smbu08_pl2frames.paf.gz')),
    ('smbu06_pl2frames', os.path.join(ATTJ, '04_controls', 'smbu06_pl2frames.paf.gz')),
]

rows = []
for ds, paf in DATASETS:
    if not os.path.exists(paf):
        print('missing', paf); continue
    stats = {k: dict(n=0, loose=0, strict=0, reads=set(), reads_strict=set()) for k in JOINTS}
    with gzip.open(paf, 'rt', encoding='utf-8', errors='replace') as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            ref = p[5]
            mapq = int(p[11])
            for jname, (jref, jpos) in JOINTS.items():
                if ref != jref:
                    continue
                ts = int(p[7])
                cg = None
                for t in p[12:]:
                    if t.startswith('cg:Z:'):
                        cg = t[5:]
                if cg is None:
                    continue
                ops = [(int(a), b) for a, b in CIG_RE.findall(cg)]
                a = anchors(ops, ts, jpos)
                if a is None:
                    continue
                stats[jname]['n'] += 1
                if mapq >= 20 and a[0] >= 1000 and a[1] >= 1000:
                    stats[jname]['loose'] += 1
                    stats[jname]['reads'].add(p[0])
                    if a[0] >= 3000 and a[1] >= 3000:
                        stats[jname]['strict'] += 1
                        stats[jname]['reads_strict'].add(p[0])
    for jname in JOINTS:
        s = stats[jname]
        rows.append([ds, jname, s['n'], s['loose'], len(s['reads']), s['strict'], len(s['reads_strict'])])
    print(ds, 'done', flush=True)

with open(os.path.join(ATTJ, '04_controls', 'state_junction_counts.tsv'), 'w', encoding='utf-8') as f:
    f.write('dataset\tjunction\talignments_spanning\tloose_alignments(mapq20,1kb)\tloose_reads\tstrict_alignments(3kb)\tstrict_reads\n')
    for r_ in rows:
        f.write('\t'.join(str(x) for x in r_) + '\n')
print('state_and_junction_counts done')
for r_ in rows:
    print(r_)
