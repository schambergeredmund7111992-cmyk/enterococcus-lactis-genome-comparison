# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""figS3_candidate_alignment.py — Supplementary Fig. S3：候选肽与已知细菌素/彼此的多序列比对"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p3lib import WORK, BLASTP, MAKEDB, read_fasta, write_fasta, run
from fig_style import plt, savefig, C
import matplotlib.gridspec as gridspec

OUT = os.path.join(WORK, '09_figures')
SEA = os.path.join(WORK, '05_bacteriocins', 'search')

prot = dict((n.split()[0], s) for n, s in read_fasta(os.path.join(SEA, 'smbu06_proteins.faa')))
C641, C679, IMM = prot['smbu06GL002641'], prot['smbu06GL002679'], prot['smbu06GL002642']


def align_pair(qn, qs, sn, ss, ev=1e-3):
    qf = os.path.join(SEA, 'tmpS3q.fa')
    sf = os.path.join(SEA, 'tmpS3s.fa')
    write_fasta(qf, [(qn, qs)])
    write_fasta(sf, [(sn, ss)])
    db = os.path.join(SEA, 'tmpS3db')
    run([MAKEDB, '-in', sf, '-dbtype', 'prot', '-out', db], log=None, check=True)
    r = run([BLASTP, '-query', qf, '-db', db, '-outfmt', '6 pident length qstart qend sstart send evalue qseq sseq',
             '-evalue', str(ev), '-seg', 'no'], log=None)
    txt = (r.stdout or '').strip()
    if not txt:
        return None
    p = txt.split('\n')[0].split('\t')   # 9 列: pident length qstart qend sstart send evalue qseq sseq
    return float(p[0]), p[7], p[8], int(p[2]), int(p[3]), int(p[4]), int(p[5])


refP = read_fasta(os.path.join(WORK, '05_bacteriocins', 'refs', 'enterocinP_O30434.fasta'))[0][1]
refH = read_fasta(os.path.join(WORK, '05_bacteriocins', 'refs', 'hiracinJM79.fasta'))[0][1]

blocks = [
    ('GL002641 vs enterocin P (O30434)', align_pair('GL002641', C641, 'enterocinP', refP)),
    ('GL002679 vs hiracin JM79', align_pair('GL002679', C679, 'hiracinJM79', refH)),
    ('GL002641 vs GL002679 (within-plasmid paralogs, low-identity alignment)', align_pair('GL002641', C641, 'GL002679', C679, ev=5.0)),
]
blocks = [(t, a) for t, a in blocks if a is not None]

fig = plt.figure(figsize=(7.2, 4.2))
gs = gridspec.GridSpec(len(blocks), 1, hspace=0.85, left=0.12)
for bi, (title, (pid, qs, ss, q0, q1, s0, s1)) in enumerate(blocks):
    ax = fig.add_subplot(gs[bi, 0])
    ax.axis('off')
    ax.set_xlim(-8, max(len(qs), len(ss)) + 2)
    ax.set_ylim(-0.5, 1.6)
    ax.text(0, 1.42, '%s — %.1f%% identity (aligned %d aa)' % (title, pid, len(qs)),
            fontsize=6.2, fontweight='bold', va='bottom')
    # 坐标串
    ax.text(-0.5, 1.05, 'cand\n%4d' % q0, fontsize=4.6, ha='right', va='center', family='DejaVu Sans Mono')
    ax.text(-0.5, 0.45, 'hmlg\n%4d' % s0, fontsize=4.6, ha='right', va='center', family='DejaVu Sans Mono')
    for j, (a, b) in enumerate(zip(qs, ss)):
        col = 'black' if a == b else 'red'
        if a != '-':
            ax.text(j, 1.05, a, fontsize=5.4, color=col, family='DejaVu Sans Mono')
        if b != '-':
            ax.text(j, 0.45, b, fontsize=5.4, color=col, family='DejaVu Sans Mono')
    # 保守 YGNGV/YDNGI 区标注
    for motif in ['YGNGV', 'YDNGI']:
        k = qs.find(motif)
        if k >= 0:
            ax.add_patch(plt.Rectangle((k - 0.4, 0.1), len(motif) - 0.2, 1.25, fill=False,
                                       edgecolor=C['bacteriocin'], lw=0.8))
            ax.text(k + 2, -0.30, motif, fontsize=5.4, color=C['bacteriocin'], ha='center')

fig.text(0.5, 0.01,
         'Red = differences. GL002641 carries the class-IIa motif variant YDNGI (vs YGNGV in enterocin P); GL002679 retains YGNGV. '
         'Sequence-level comparison only — no activity data.',
         fontsize=5.8, ha='center')

savefig(fig, os.path.join(OUT, 'FigS3_candidate_alignment', 'FigS3'))
print('figS3 done')
