# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
junction_support2.py — 03_supporting_reads（本会话独立实现）
对 panel_db（13 条参考）与 full_db（4 条复制子参考）的 SAM/PAF 做双路径解析：

预注册判据（分析前固定，见代码常量 CRIT）：
  单条 read 的一条 alignment：
    - junction 位置被 M/=/X（匹配）覆盖（不可是删除）
    - 两侧锚定 = junction 两侧在比对内的匹配参考碱基数
    - loose: 两侧 >=1000 bp；strict: 两侧 >=3000 bp
    - mapQ >= 20；primary alignment（非 secondary/supplementary）
输出：
  03_supporting_reads/junction_candidates2.tsv       panel 全部 junction 相关 alignment（含未过判据）
  03_supporting_reads/junction_support_summary.tsv    各 junction 候选/通过数（SAM 与 PAF 两路径）
  03_supporting_reads/alignment_evidence_<ref>.tsv    通过 loose 的每条 alignment 证据
  03_supporting_reads/supporting_reads_<ref>.fasta    通过 strict 的 read 序列（原始方向）
  03_supporting_reads/read_best_by_ref.tsv            full_db: 每 read 对每参考的最佳 AS/mapQ（竞争参考分析用）
"""
import gzip
import os
import re
from collections import defaultdict

ATTJ = r'<PROJECT_ROOT>\longread_validation'
MAP = os.path.join(ATTJ, '02_mapping')
SR = os.path.join(ATTJ, '03_supporting_reads')
FQ = r'<PROJECT_ROOT>\smbu08\smbu08\1.Cleandata\smbu08.filtered_reads.fq.gz'
LOG = open(os.path.join(ATTJ, '00_manifest', 'logs', 'junction_support2.log'), 'a', encoding='utf-8')

CRIT = dict(min_mapq=20, anchor_loose=1000, anchor_strict=3000)

# junction 定义
jdef = {}
for line in open(os.path.join(ATTJ, '06_tables', 'junction_definitions.tsv'), encoding='utf-8').read().splitlines()[1:]:
    p = line.split('\t')
    jdef[p[0]] = (p[1], int(p[2]))
LOG.write('junction definitions: %s\n' % jdef)

CIG_RE = re.compile(r'(\d+)([MIDNSHP=X])')


def parse_cigar(c):
    return [(int(n), o) for n, o in CIG_RE.findall(c)]


def anchor_info(ops, r0, jpos):
    """返回 dict(covered, left, right, maxD_left, maxD_right, maxI_total) 或 None（junction 不在匹配内）"""
    pos = r0
    left = right = 0
    maxD_l = maxD_r = 0
    ins_tot = 0
    found = False
    for n, o in ops:
        if o in 'M=X':          # 注意：minimap2 默认 CIGAR 使用 M（非 =/X）
            if pos <= jpos < pos + n:
                found = True
                left += jpos - pos
                right += pos + n - jpos - 1
            elif pos + n <= jpos:
                left += n
            else:
                right += n
            pos += n
        elif o == 'I':
            ins_tot += n
        elif o in 'DN':
            if pos <= jpos < pos + n:
                return None
            if pos + n <= jpos:
                maxD_l = max(maxD_l, n)
            else:
                maxD_r = max(maxD_r, n)
            pos += n
        elif o in 'SH':
            pass
    if not found:
        return None
    return dict(left=left, right=right, maxD_l=maxD_l, maxD_r=maxD_r, ins=ins_tot)


def ident_from(cigar_ops, nm, aln_len_ref):
    m = sum(n for n, o in cigar_ops if o in 'M=X')
    ins = sum(n for n, o in cigar_ops if o == 'I')
    dele = sum(n for n, o in cigar_ops if o in 'DN')
    mism = nm - ins - dele if nm is not None and nm >= 0 else None
    if mism is None or m == 0:
        return None
    return (m - mism) / m


# ---------------- 1) panel SAM 解析 ----------------
cand_rows = []
summary = defaultdict(lambda: defaultdict(int))
support_strict = defaultdict(set)
support_loose = defaultdict(set)

with gzip.open(os.path.join(MAP, 'panel_db.sam.gz'), 'rt', encoding='utf-8', errors='replace') as f:
    for line in f:
        if line[0] == '@':
            continue
        p = line.rstrip('\n').split('\t')
        flag = int(p[1]); ref = p[2]
        if ref == '*':
            continue
        summary[ref]['alignments'] += 1
        is_sec = bool(flag & 0x100); is_sup = bool(flag & 0x800)
        primary = not is_sec and not is_sup
        if primary:
            summary[ref]['primary'] += 1
        mapq = int(p[4])
        ops = parse_cigar(p[5])
        r0 = int(p[3]) - 1
        tags = {}
        for t in p[11:]:
            k, _, v = t.split(':', 2)
            tags[k] = v
        nm = int(tags['NM']) if 'NM' in tags else None
        asc = int(tags['AS']) if 'AS' in tags else None
        rlen = len(p[9]) if p[9] != '*' else 0
        if ref not in jdef:
            continue
        jtype, jpos = jdef[ref]
        a = anchor_info(ops, r0, jpos)
        covered = a is not None
        if covered:
            summary[ref]['spanning'] += 1
            if primary:
                summary[ref]['spanning_primary'] += 1
            if primary and mapq >= CRIT['min_mapq']:
                summary[ref]['spanning_primary_mapq20'] += 1
        lo = bool(a and a['left'] >= CRIT['anchor_loose'] and a['right'] >= CRIT['anchor_loose'])
        st = bool(a and a['left'] >= CRIT['anchor_strict'] and a['right'] >= CRIT['anchor_strict'])
        ok_mq = mapq >= CRIT['min_mapq']
        ident = ident_from(ops, nm, None)
        clip_l = ops[0][0] if ops and ops[0][1] == 'S' else 0
        clip_r = ops[-1][0] if ops and ops[-1][1] == 'S' else 0
        rs = sum(n for n, o in ops if o in 'M=XDN')
        if primary and lo and ok_mq:
            support_loose[ref].add(p[0])
        if primary and st and ok_mq:
            support_strict[ref].add(p[0])
        cand_rows.append(dict(read=p[0], ref=ref, jtype=jtype, flag=flag, mapq=mapq, primary=primary,
                              AS=asc, NM=nm, ident=round(ident, 4) if ident is not None else '',
                              r0=r0, r1=r0 + rs, clip_l=clip_l, clip_r=clip_r, read_len=rlen,
                              a_left=a['left'] if a else '', a_right=a['right'] if a else '',
                              maxD_l=a['maxD_l'] if a else '', maxD_r=a['maxD_r'] if a else '',
                              ins=a['ins'] if a else '', covered=covered,
                              pass_loose=bool(lo and primary and ok_mq),
                              pass_strict=bool(st and primary and ok_mq), cigar=p[5]))

cols = ['read', 'ref', 'jtype', 'flag', 'mapq', 'primary', 'AS', 'NM', 'ident', 'r0', 'r1',
        'clip_l', 'clip_r', 'read_len', 'a_left', 'a_right', 'maxD_l', 'maxD_r', 'ins',
        'covered', 'pass_loose', 'pass_strict', 'cigar']
with open(os.path.join(SR, 'junction_candidates2.tsv'), 'w', encoding='utf-8') as f:
    f.write('\t'.join(cols) + '\n')
    for c in sorted(cand_rows, key=lambda x: (x['ref'], not x['pass_strict'], -x['mapq'])):
        f.write('\t'.join(str(c[k]) for k in cols) + '\n')

# ---------------- 2) PAF 交叉核验（独立格式/解析路径） ----------------
paf_summary = defaultdict(lambda: defaultdict(int))
paf_support = defaultdict(set)
with gzip.open(os.path.join(MAP, 'panel_db.paf.gz'), 'rt', encoding='utf-8', errors='replace') as f:
    for line in f:
        p = line.rstrip('\n').split('\t')
        ref = p[5]
        if ref not in jdef:
            continue
        qname, qlen, qs, qe, strand = p[0], int(p[1]), int(p[2]), int(p[3]), p[4]
        tlen, ts, te = int(p[6]), int(p[7]), int(p[8])
        mapq = int(p[11])
        cg = None
        for t in p[12:]:
            if t.startswith('cg:Z:'):
                cg = t[5:]
        if cg is None:
            continue
        ops = parse_cigar(cg)
        jtype, jpos = jdef[ref]
        paf_summary[ref]['alignments'] += 1
        a = anchor_info(ops, ts, jpos)
        if a:
            paf_summary[ref]['spanning'] += 1
            if mapq >= CRIT['min_mapq']:
                paf_summary[ref]['spanning_mapq20'] += 1
                if a['left'] >= CRIT['anchor_loose'] and a['right'] >= CRIT['anchor_loose']:
                    paf_summary[ref]['loose_mapq20'] += 1
                if a['left'] >= CRIT['anchor_strict'] and a['right'] >= CRIT['anchor_strict']:
                    paf_summary[ref]['strict_mapq20'] += 1
                    paf_support[ref].add(qname)

with open(os.path.join(SR, 'junction_support_summary.tsv'), 'w', encoding='utf-8') as f:
    f.write('ref\tjtype\tjpos\tSAM_alignments\tSAM_spanning\tSAM_spanning_primary\tSAM_spanning_primary_mapq20\t'
            'SAM_loose_support_reads\tSAM_strict_support_reads\tPAF_alignments\tPAF_spanning\tPAF_spanning_mapq20\t'
            'PAF_loose_mapq20\tPAF_strict_mapq20\n')
    for ref, (jtype, jpos) in jdef.items():
        s = summary[ref]; ps = paf_summary[ref]
        f.write('%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n' % (
            ref, jtype, jpos, s['alignments'], s['spanning'], s['spanning_primary'],
            s['spanning_primary_mapq20'], len(support_loose[ref]), len(support_strict[ref]),
            ps['alignments'], ps['spanning'], ps['spanning_mapq20'], ps['loose_mapq20'], ps['strict_mapq20']))
LOG.write('SAM vs PAF summary written\n')

# 一致性检查：两路径 strict 支持 read 集合差异
for ref in jdef:
    only_sam = support_strict[ref] - paf_support[ref] - {''}
    only_paf = paf_support[ref] - support_strict[ref] - {''}
    LOG.write('strict support %s: SAM=%d PAF=%d  SAM_only=%d PAF_only=%d\n'
              % (ref, len(support_strict[ref]), len(paf_support[ref]), len(only_sam), len(only_paf)))

# ---------------- 3) 证据表 + 支持 reads FASTA ----------------
wanted = set()
for ref in ['attJ_junction', 'unit1_to_unit2_junction', 'unit2_to_unit1_circular_junction',
            'retained_window'] + [r for r in jdef if r.startswith('neg_')]:
    wanted |= support_loose.get(ref, set())

seqs = {}
if wanted:
    with gzip.open(FQ, 'rt') as f:
        while True:
            h = f.readline()
            if not h:
                break
            s = f.readline().strip()
            f.readline(); f.readline()
            nm = h[1:].split()[0]
            if nm in wanted:
                seqs[nm] = s
print('candidate reads extracted:', len(seqs), flush=True)

for ref in jdef:
    ev = [c for c in cand_rows if c['ref'] == ref and c['pass_loose']]
    if ev:
        ev.sort(key=lambda x: (-(x['a_left'] + x['a_right'])))
        with open(os.path.join(SR, 'alignment_evidence_%s.tsv' % ref), 'w', encoding='utf-8') as f:
            f.write('\t'.join(cols) + '\n')
            for c in ev:
                f.write('\t'.join(str(c[k]) for k in cols) + '\n')
    for tag, sset in [('strict', support_strict[ref]), ('loose', support_loose[ref])]:
        if sset:
            with open(os.path.join(SR, 'supporting_reads_%s_%s.fasta' % (ref, tag)), 'w') as f:
                for n in sorted(sset):
                    if n in seqs:
                        f.write('>%s\n%s\n' % (n, seqs[n]))

# ---------------- 4) full_db：read × ref 最佳 AS（竞争参考分析输入） ----------------
best = defaultdict(dict)
with gzip.open(os.path.join(MAP, 'full_db.sam.gz'), 'rt', encoding='utf-8', errors='replace') as f:
    for line in f:
        if line[0] == '@':
            continue
        p = line.rstrip('\n').split('\t')
        ref = p[2]
        if ref == '*':
            continue
        flag = int(p[1])
        asc = None; nm = None
        for t in p[11:]:
            if t.startswith('AS:i:'):
                asc = int(t[5:])
            elif t.startswith('NM:i:'):
                nm = int(t[5:])
        if asc is None:
            continue
        ops = parse_cigar(p[5])
        rs = sum(n for n, o in ops if o in 'M=XDN')
        prim = not (flag & 0x100) and not (flag & 0x800)
        rec = best[p[0]]
        cur = rec.get(ref)
        if cur is None or asc > cur[0]:
            rec[ref] = (asc, int(p[4]), 1 if prim else 0, int(p[3]) - 1, int(p[3]) - 1 + rs, flag)
with open(os.path.join(SR, 'read_best_by_ref.tsv'), 'w', encoding='utf-8') as f:
    f.write('read\tref\tAS\tmapq\tprimary\tr0\tr1\tflag\n')
    for rd, refs in best.items():
        for ref, v in refs.items():
            f.write('%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n' % (rd, ref, v[0], v[1], v[2], v[3], v[4], v[5]))
LOG.write('read_best_by_ref rows=%d reads=%d\n' % (sum(len(v) for v in best.values()), len(best)))

print('junction_support2 done; candidates=%d' % len(cand_rows))
