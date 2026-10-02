# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
physchem.py — 候选细菌素成熟肽理化性质与同源对齐（模块 05 补充）
- 与已知细菌素（enterocin P O30434；hiracin JM79）做对齐级比较（blastp qseq/sseq）
- 理化性质：长度/MW/pI/GRAVY/净电荷（简单 pKa 模型，记录为计算值）
输出：05_bacteriocins/search/physchem_candidates.tsv, alignment_*.txt
"""
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
from p3lib import WORK, BLASTP, MAKEDB, read_fasta, write_fasta, run, open_log

OUT = os.path.join(WORK, '05_bacteriocins')
SEA = os.path.join(OUT, 'search')
LOG = open_log(os.path.join(OUT, 'logs', 'physchem.log'))

# 候选序列（来自 GBK 翻译，与 candidate_screen 相同来源）
prot = dict((n.split()[0], s) for n, s in read_fasta(os.path.join(SEA, 'smbu06_proteins.faa')))
C641 = prot['smbu06GL002641']
C679 = prot['smbu06GL002679']
IMM = prot['smbu06GL002642']
LOG.write('641= %s\n679= %s\n642= %s\n' % (C641, C679, IMM))

# 已知参考：fetch O30434 (enterocin P) 与 hiracin JM79（UniProt 搜索）
REF = os.path.join(OUT, 'refs')
refs = {}
for name, url in [('enterocinP_O30434', 'https://rest.uniprot.org/uniprotkb/O30434.fasta'),
                  ('hiracinJM79', 'https://rest.uniprot.org/uniprotkb/search?query=hiracin+JM79&format=fasta&size=3')]:
    r = run(['curl', '-sL', '--max-time', '60', url], log=LOG)
    if r.stdout.startswith('>'):
        p = os.path.join(REF, name + '.fasta')
        with open(p, 'w') as f:
            f.write(r.stdout)
        recs = read_fasta(p)
        refs[name] = recs[0]
        LOG.write('ref %s: %s len=%d\n' % (name, recs[0][0][:60], len(recs[0][1])))

# 对齐（blastp qseq/sseq）
def align_pair(qname, qseq, sname, sseq):
    qf = os.path.join(SEA, 'tmp_q.fa')
    sf = os.path.join(SEA, 'tmp_s.fa')
    write_fasta(qf, [(qname, qseq)])
    write_fasta(sf, [(sname, sseq)])
    db = os.path.join(SEA, 'tmp_s_db')
    run([MAKEDB, '-in', sf, '-dbtype', 'prot', '-out', db], log=None, check=True)
    r = run([BLASTP, '-query', qf, '-db', db, '-outfmt', '6 pident length qstart qend sstart send evalue qseq sseq',
             '-evalue', '1e-3', '-seg', 'no'], log=None)
    return r.stdout.strip()


AA_MW = {'A': 71.08, 'R': 156.19, 'N': 114.10, 'D': 115.09, 'C': 103.14, 'E': 129.12, 'Q': 128.13,
         'G': 57.05, 'H': 137.14, 'I': 113.16, 'L': 113.16, 'K': 128.17, 'M': 131.19, 'F': 147.18,
         'P': 97.12, 'S': 87.08, 'T': 101.10, 'W': 186.21, 'Y': 163.18, 'V': 99.13}
KD = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'E': -3.5, 'Q': -3.5, 'G': -0.4, 'H': -3.2,
      'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9,
      'Y': -1.3, 'V': 4.2}
PKA = {'Cterm': 3.55, 'Nterm': 7.5, 'D': 4.05, 'E': 4.45, 'C': 9.0, 'Y': 10.0, 'H': 5.98, 'K': 10.0, 'R': 12.0}


def mw(s):
    return sum(AA_MW.get(a, 0) for a in s) + 18.02


def gravy(s):
    return sum(KD.get(a, 0) for a in s) / len(s)


def charge(s, pH):
    q = 0.0
    q += -1 / (1 + 10 ** (PKA['Cterm'] - pH))
    q += 1 / (1 + 10 ** (pH - PKA['Nterm']))
    for a, pka, sign in [('D', PKA['D'], -1), ('E', PKA['E'], -1), ('C', PKA['C'], -1),
                         ('Y', PKA['Y'], -1), ('H', PKA['H'], 1), ('K', PKA['K'], 1), ('R', PKA['R'], 1)]:
        n = s.count(a)
        if n:
            if sign < 0:
                q += -n / (1 + 10 ** (pka - pH))
            else:
                q += n / (1 + 10 ** (pH - pka))
    return q


def pI(s):
    lo, hi = 2.0, 13.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if charge(s, mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


rows = []
defs = [('GL002641 full (enterocin P-like)', C641),
        ('GL002641 last-40aa (predicted mature-like region by homolog alignment)', C641[-40:]),
        ('GL002679 full (hiracin JM79-like)', C679),
        ('GL002679 last-44aa (mature-like region)', C679[-44:]),
        ('GL002642 full (Sakacin-A immunity factor)', IMM)]
for name, s in defs:
    rows.append(dict(seq_name=name, length_aa=len(s), mw_da=round(mw(s), 1), pI=round(pI(s), 2),
                     gravy=round(gravy(s), 3), charge_pH7=round(charge(s, 7.0), 2),
                     n_cys=s.count('C'), sequence=s))

with open(os.path.join(SEA, 'physchem_candidates.tsv'), 'w', encoding='utf-8') as f:
    f.write('seq_name\tlength_aa\tmw_da\tpI\tgravy\tcharge_pH7\tn_cys\tsequence\n')
    for r in rows:
        f.write('\t'.join(str(r[k]) for k in ['seq_name', 'length_aa', 'mw_da', 'pI', 'gravy',
                                              'charge_pH7', 'n_cys', 'sequence']) + '\n')
    LOG.write('%s\n' % r)

# 对齐输出
if 'enterocinP_O30434' in refs:
    aln = align_pair('GL002641', C641, 'enterocinP_O30434', refs['enterocinP_O30434'][1])
    with open(os.path.join(SEA, 'alignment_GL002641_vs_enterocinP.txt'), 'w') as f:
        f.write('# GL002641 vs enterocin P (O30434), blastp -seg no\n#' + aln + '\n')
    LOG.write('align 641 vs enterocin P: %s\n' % aln.split('\t')[:8])
if 'hiracinJM79' in refs:
    aln2 = align_pair('GL002679', C679, 'hiracinJM79', refs['hiracinJM79'][1])
    with open(os.path.join(SEA, 'alignment_GL002679_vs_hiracinJM79.txt'), 'w') as f:
        f.write('# GL002679 vs hiracin JM79, blastp -seg no\n#' + aln2 + '\n')
    LOG.write('align 679 vs hiracin: %s\n' % aln2.split('\t')[:8])
print('physchem done')
