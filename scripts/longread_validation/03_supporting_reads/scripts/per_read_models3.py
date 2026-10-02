# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
per_read_models3.py — 扩展候选集（loose 支持 ∪ 二聚体 junction primary 跨越 reads）的
每 read × 每模型竞争比对与 strong/ambiguous/weak 判定（预注册阈值同 models2）。
输出：03_supporting_reads/read_classification3.tsv, model_scores3.tsv
"""
import gzip
import os
import subprocess
from collections import defaultdict

ATTJ = r'<PROJECT_ROOT>\longread_validation'
REFDIR = os.path.join(ATTJ, '01_references', 'references')
MM2 = r'<TOOLS_ROOT>\minimap2-win\minimap2-2.31-r1302-windows-x86_64-ucrt64\minimap2.exe'
SR = os.path.join(ATTJ, '03_supporting_reads')
FQ = r'<PROJECT_ROOT>\smbu08\smbu08\1.Cleandata\smbu08.filtered_reads.fq.gz'
LOG = open(os.path.join(ATTJ, '00_manifest', 'logs', 'per_read_models3.log'), 'a', encoding='utf-8')

MODELS = [
    ('attJ_junction', 'hypothesis_attJ'),
    ('chr08_full', 'hypothesis_attJ'),
    ('retained_window', 'alternative_retained'),
    ('chr06_full', 'alternative_retained'),
    ('unit1_to_unit2_junction', 'hypothesis_dimer_fwd'),
    ('unit2_to_unit1_circular_junction', 'hypothesis_dimer_back'),
    ('pl2_asassembled', 'model_linear_dimer'),
    ('dimer_back_internal', 'model_circular_dimer'),
    ('monomer_circle_repA', 'model_monomer'),
    ('monomer_circle_repB', 'model_monomer'),
]

# 收集候选（读 junction_candidates2.tsv）
hdr = None
cand_reads = {}      # read -> set(jtype)
with open(os.path.join(SR, 'junction_candidates2.tsv'), encoding='utf-8') as f:
    hdr = f.readline().rstrip('\n').split('\t')
    for line in f:
        p = dict(zip(hdr, line.rstrip('\n').split('\t')))
        take = False
        if p['pass_loose'] == 'True':
            take = True
        if p['ref'] in ('unit1_to_unit2_junction', 'unit2_to_unit1_circular_junction') \
           and p['primary'] == 'True' and p['covered'] == 'True':
            take = True
        if take:
            cand_reads.setdefault(p['read'], set()).add(p['jtype'])
LOG.write('candidate reads: %d\n' % len(cand_reads))

# 序列：先取已有 loose FASTA，再从 fq.gz 补齐
seqs = {}
for fn in os.listdir(SR):
    if fn.startswith('supporting_reads_') and fn.endswith('_loose.fasta'):
        name = None
        for line in open(os.path.join(SR, fn)):
            if line[0] == '>':
                name = line[1:].strip()
            else:
                seqs[name] = line.strip()
missing = set(cand_reads) - set(seqs)
if missing:
    with gzip.open(FQ, 'rt') as f:
        while True:
            h = f.readline()
            if not h:
                break
            s = f.readline().strip()
            f.readline(); f.readline()
            nm = h[1:].split()[0]
            if nm in missing:
                seqs[nm] = s
LOG.write('seqs: %d (补齐 %d)\n' % (len(seqs), len(missing)))
cand_fa = os.path.join(SR, 'candidate_reads_expanded.fasta')
with open(cand_fa, 'w') as f:
    for n, s in seqs.items():
        f.write('>%s\n%s\n' % (n, s))

scores = defaultdict(dict)
for ref, mtype in MODELS:
    rf = os.path.join(REFDIR, ref + '.fasta')
    r = subprocess.run([MM2, '-c', '-x', 'map-ont', '-t', '12', '--secondary=no', rf, cand_fa],
                       capture_output=True, text=True, encoding='utf-8', errors='replace')
    best = {}
    for line in r.stdout.strip().split('\n'):
        if not line:
            continue
        p = line.split('\t')
        q = p[0]
        mlen = int(p[9]); alen = int(p[10])
        AS = 0
        for t in p[12:]:
            if t.startswith('AS:i:'):
                AS = int(t[5:])
        ident = mlen / alen if alen else 0
        if q not in best or AS > best[q][0]:
            best[q] = (AS, ident, alen, p[4])
    for q, (AS, ident, alen, strand) in best.items():
        scores[q][ref] = dict(AS=AS, ident=round(ident, 4), alen=alen, strand=strand)
    print('%s: %d/%d' % (ref, len(best), len(seqs)), flush=True)

rows = []
for read, models in scores.items():
    ranked = sorted(models.items(), key=lambda kv: -kv[1]['AS'])
    best_ref, best_v = ranked[0]
    second_ref, second_v = (ranked[1] if len(ranked) > 1 else (None, dict(AS=0, ident=0, alen=0)))
    delta = best_v['AS'] - second_v['AS']
    rlen = len(seqs.get(read, ''))
    cover = best_v['alen'] / rlen if rlen else 0
    strong = (delta >= 10) and (best_v['ident'] >= 0.90) and (cover >= 0.30)
    ambiguous = (delta < 10) or (best_v['ident'] < 0.90)
    verdict = 'strong' if strong else ('ambiguous' if ambiguous else 'weak')
    rows.append([read, ','.join(sorted(cand_reads.get(read, {'?'}))), rlen, best_ref, best_v['AS'],
                 second_ref, second_v['AS'], delta, best_v['ident'], second_v['ident'],
                 round(cover, 3), verdict])

with open(os.path.join(SR, 'read_classification3.tsv'), 'w', encoding='utf-8') as f:
    f.write('read\tjtype\tread_len\tbest_model\tbest_AS\tsecond_model\tsecond_AS\tdelta\tbest_ident\tsecond_ident\tcover\tverdict\n')
    for r_ in sorted(rows, key=lambda x: (-x[4])):
        f.write('\t'.join(str(x) for x in r_) + '\n')
with open(os.path.join(SR, 'model_scores3.tsv'), 'w', encoding='utf-8') as f:
    f.write('read\tmodel\tAS\tidentity\taln_len\tstrand\n')
    for read, models in sorted(scores.items()):
        for m, v in sorted(models.items(), key=lambda kv: -kv[1]['AS']):
            f.write('%s\t%s\t%d\t%s\t%d\t%s\n' % (read, m, v['AS'], v['ident'], v['alen'], v['strand']))

from collections import Counter
LOG.write('verdicts: %s\n' % dict(Counter(r_[11] for r_ in rows)))
LOG.write('best_model: %s\n' % dict(Counter(r_[3] for r_ in rows)))
print('n=%d verdicts:' % len(rows), dict(Counter(r_[11] for r_ in rows)))
print('best_model:', dict(Counter(r_[3] for r_ in rows)))
