# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
per_read_models2.py — 每条候选 read 对各结构模型单独做一次比对（minimap2 -c，--secondary=no），
取每模型最佳 AS/identity/覆盖长度，比较最优与次优（预注册阈值）。
模型集合（每个模型一个参考序列，避免同一参考多拷贝造成的内部竞争）：
  hypothesis_attJ      : attJ_junction, chr08_full
  alternative_retained : retained_window, chr06_full
  hypothesis_dimer_fwd : unit1_to_unit2_junction
  hypothesis_dimer_back: unit2_to_unit1_circular_junction
  model_linear_dimer   : pl2_asassembled
  model_circular_dimer : dimer_back_internal
  model_monomer        : monomer_circle_repA, monomer_circle_repB
预注册判定：
  strong    : best AS - second AS >= 10 且 identity >= 0.90 且 aln 覆盖 >= 30% read
  ambiguous : delta < 10 或 identity < 0.90
  weak      : 其余（覆盖不足等）
输出：
  03_supporting_reads/model_scores2.tsv
  03_supporting_reads/read_classification2.tsv
"""
import gzip
import os
import subprocess
from collections import defaultdict

ATTJ = r'<PROJECT_ROOT>\longread_validation'
REFDIR = os.path.join(ATTJ, '01_references', 'references')
MM2 = r'<TOOLS_ROOT>\minimap2-win\minimap2-2.31-r1302-windows-x86_64-ucrt64\minimap2.exe'
SR = os.path.join(ATTJ, '03_supporting_reads')
LOG = open(os.path.join(ATTJ, '00_manifest', 'logs', 'per_read_models2.log'), 'a', encoding='utf-8')

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

# ---- 1) 候选 reads：合并所有 loose 支持 FASTA ----
seqs = {}
for fn in os.listdir(SR):
    if fn.startswith('supporting_reads_') and fn.endswith('_loose.fasta'):
        name = None
        for line in open(os.path.join(SR, fn)):
            if line[0] == '>':
                name = line[1:].strip()
            else:
                seqs[name] = line.strip()
LOG.write('loose candidate reads (union): %d\n' % len(seqs))
cand_fa = os.path.join(SR, 'candidate_reads_loose.fasta')
with open(cand_fa, 'w') as f:
    for n, s in seqs.items():
        f.write('>%s\n%s\n' % (n, s))

# ---- 2) 每模型一次比对 ----
scores = defaultdict(dict)
for ref, mtype in MODELS:
    rf = os.path.join(REFDIR, ref + '.fasta')
    cmd = [MM2, '-c', '-x', 'map-ont', '-t', '12', '--secondary=no', rf, cand_fa]
    r = subprocess.run(cmd, capture_output=True, text=True, encoding='utf-8', errors='replace')
    LOG.write('model %s rc=%d\n' % (ref, r.returncode))
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
    print('%s: %d/%d mapped' % (ref, len(best), len(seqs)), flush=True)

# ---- 3) 判定 ----
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
    rows.append([read, rlen, best_ref, best_v['AS'], second_ref, second_v['AS'], delta,
                 best_v['ident'], second_v['ident'], round(cover, 3), verdict])

with open(os.path.join(SR, 'read_classification2.tsv'), 'w', encoding='utf-8') as f:
    f.write('read\tread_len\tbest_model\tbest_AS\tsecond_model\tsecond_AS\tdelta\tbest_ident\tsecond_ident\tcover\tverdict\n')
    for r_ in sorted(rows, key=lambda x: (-x[3])):
        f.write('\t'.join(str(x) for x in r_) + '\n')

with open(os.path.join(SR, 'model_scores2.tsv'), 'w', encoding='utf-8') as f:
    f.write('read\tmodel\tAS\tidentity\taln_len\tstrand\n')
    for read, models in sorted(scores.items()):
        for m, v in sorted(models.items(), key=lambda kv: -kv[1]['AS']):
            f.write('%s\t%s\t%d\t%s\t%d\t%s\n' % (read, m, v['AS'], v['ident'], v['alen'], v['strand']))

from collections import Counter
LOG.write('verdicts: %s\n' % dict(Counter(r_[10] for r_ in rows)))
LOG.write('best_model counts: %s\n' % dict(Counter(r_[2] for r_ in rows)))
print('per_read_models2 done; n=%d' % len(rows))
print('best_model counts:', dict(Counter(r_[2] for r_ in rows)))
