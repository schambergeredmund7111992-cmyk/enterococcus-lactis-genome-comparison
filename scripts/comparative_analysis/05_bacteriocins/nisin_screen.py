# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
nisin_screen.py — 模块 05：nisin 双层检索（蛋白 tblastn + 核酸 blastn）+ 阳性/阴性对照
目标：smbu06 完整组装、smbu08 完整组装、F44 阳性对照簇区段、IL1403 阴性对照基因组
判据（沿用已验证方法）：完整簇 = nisA+nisB+nisC 均命中且 identity≥80% 且 qcov≥80%
输出：05_bacteriocins/search/nisin_tblastn_hits.tsv, nisin_status_matrix.tsv, nisin_control_validation.tsv
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
from p3lib import (WORK, IN06, IN08, BLASTN, BLASTP, TBLASTN, MAKEDB, read_fasta,
                   write_fasta, run, open_log)

OUT = os.path.join(WORK, '05_bacteriocins')
REF = os.path.join(OUT, 'refs')
SEA = os.path.join(OUT, 'search')
os.makedirs(SEA, exist_ok=True)
LOG = open_log(os.path.join(OUT, 'logs', 'nisin_screen.log'))

NISIN_FAA = os.path.join(REF, 'nisin_reference_proteins.faa')
F44 = os.path.join(REF, 'nisin_positive_cluster_F44.fasta')
IL1403 = WORK + r'\02_taxonomy\refs\GCF_000006865.1\ncbi_dataset\data\GCF_000006865.1\GCF_000006865.1_ASM686v1_genomic.fna'

TARGETS = {
    'smbu06_complete': IN06 + r'\smbu06\2.Assembly\smbu06.Complete.genome.fasta',
    'smbu08_complete': IN08 + r'\smbu08\2.Assembly\smbu08.Complete.genome.fasta',
    'F44_positive_control': F44,
    'IL1403_negative_control': IL1403,
}
META = ['pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send',
        'evalue', 'bitscore', 'qlen', 'slen']


def make_nt_db(path, name):
    db = os.path.join(SEA, name + '_db')
    run([MAKEDB, '-in', path, '-dbtype', 'nucl', '-out', db], log=LOG, check=True)
    return db


hits = []
for tname, tpath in TARGETS.items():
    db = make_nt_db(tpath, tname)
    for tool, query, qkind, extra in [(TBLASTN, NISIN_FAA, 'prot', ['-seg', 'no']),
                                      (BLASTN, F44, 'nucl', [])]:
        cmd = [tool, '-query', query, '-db', db,
               '-outfmt', '6 qseqid sseqid ' + ' '.join(META),
               '-evalue', '1e-5', '-max_target_seqs', '5', '-max_hsps', '50', '-num_threads', '8'] + extra
        r = run(cmd, log=LOG)
        if r.returncode != 0:
            LOG.write('!! %s vs %s rc=%d\n' % (tool, tname, r.returncode))
            continue
        for line in r.stdout.strip().split('\n'):
            if not line:
                continue
            p = line.split('\t')
            qid, sid = p[0], p[1]
            vals = dict(zip(META, p[2:]))
            qcov = (int(vals['qend']) - int(vals['qstart']) + 1) / float(vals['qlen']) * 100
            hits.append(dict(search=('tblastn' if tool == TBLASTN else 'blastn'), target=tname,
                             query=qid, subject=sid,
                             pident=float(vals['pident']), alen=int(vals['length']),
                             evalue=float(vals['evalue']), qcov=round(qcov, 1),
                             sstart=vals['sstart'], send=vals['send']))
        LOG.write('%s (%s) vs %s: %d hits\n' % (tool.split('\\')[-1], qkind, tname,
                                                sum(1 for h in hits if h['target'] == tname and h['search'] == ('tblastn' if tool == TBLASTN else 'blastn'))))

with open(os.path.join(SEA, 'nisin_tblastn_hits.tsv'), 'w', encoding='utf-8') as f:
    f.write('search\ttarget\tquery\tsubject\tpident\talen\tevalue\tqcov_pct\tsstart\tsend\n')
    for h in hits:
        f.write('\t'.join(str(h[k]) for k in ['search', 'target', 'query', 'subject', 'pident', 'alen',
                                              'evalue', 'qcov', 'sstart', 'send']) + '\n')

# 状态矩阵（tblastn，逐蛋白最佳命中）
NIS_SPECIFIC = ['nisA', 'nisB', 'nisC', 'nisI', 'nisP', 'nisR']
rows = []
for strain in ['smbu06_complete', 'smbu08_complete', 'F44_positive_control', 'IL1403_negative_control']:
    sh = [h for h in hits if h['target'] == strain and h['search'] == 'tblastn']
    best = {}
    for h in sh:
        q = h['query'].split('|')[0]
        if q not in best or (h['pident'], h['qcov']) > (best[q]['pident'], best[q]['qcov']):
            best[q] = h
    specific_hits = [q for q in NIS_SPECIFIC if q in best and best[q]['pident'] >= 80 and best[q]['qcov'] >= 80]
    status = 'complete_cluster' if all(q in specific_hits for q in ['nisA', 'nisB', 'nisC']) else \
             ('partial' if specific_hits else 'not_detected')
    rows.append(dict(strain=strain, status=status,
                     details=';'.join('%s:%s%%/%s%%' % (q, best[q]['pident'], best[q]['qcov'])
                                      for q in sorted(best))))
    LOG.write('STATUS %s: %s | %s\n' % (strain, status, rows[-1]['details']))

with open(os.path.join(SEA, 'nisin_status_matrix.tsv'), 'w', encoding='utf-8') as f:
    f.write('strain\tstatus\tdetails(nisin蛋白最佳命中 identity/qcov)\n')
    for r in rows:
        f.write('%s\t%s\t%s\n' % (r['strain'], r['status'], r['details']))

# 对照验证表
ctrl = [r for r in rows if 'control' in r['strain']]
ok_pos = any(r['strain'] == 'F44_positive_control' and r['status'] == 'complete_cluster' for r in ctrl)
ok_neg = any(r['strain'] == 'IL1403_negative_control' and r['status'] == 'not_detected' for r in ctrl)
with open(os.path.join(SEA, 'nisin_control_validation.tsv'), 'w', encoding='utf-8') as f:
    f.write('control\texpectation\tobserved\tpass\n')
    f.write('F44_positive_control\tcomplete nisin cluster\t%s\t%s\n' % (
        [r['status'] for r in ctrl if r['strain'] == 'F44_positive_control'][0], ok_pos))
    f.write('IL1403_negative_control\tno nisin\t%s\t%s\n' % (
        [r['status'] for r in ctrl if r['strain'] == 'IL1403_negative_control'][0], ok_neg))
print('nisin_screen done; control validation:', ok_pos, ok_neg)
