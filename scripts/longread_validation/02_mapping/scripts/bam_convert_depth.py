# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
bam_convert_depth.py — 用 samtools（MSYS2 1.24）把 SAM.gz 转 coordinate-sorted BAM + 索引，
  flagstat/idxstats QC，并对关键区域做 samtools depth（深度算法A，独立实现），
  同时输出 junction ±2kb 的 BAM 子集（供图表/复核）。
输出：02_mapping/{full,panel}_db.sorted.bam(+.bai), *.flagstat.txt, *.idxstats.txt
      02_mapping/bam_subsets/*.bam
      04_controls/depth_samtools_regions.tsv
"""
import os
import subprocess

ATTJ = r'<PROJECT_ROOT>\longread_validation'
MAP = os.path.join(ATTJ, '02_mapping')
OUT = os.path.join(ATTJ, '04_controls')
ST = os.path.join(ATTJ, '00_manifest', 'tools', 'bin', 'samtools.exe')
LOG = open(os.path.join(ATTJ, '00_manifest', 'logs', 'bam_convert.log'), 'a', encoding='utf-8')

V = {}
for line in open(os.path.join(ATTJ, '01_references', 'independent_verification.tsv'), encoding='utf-8').read().splitlines()[1:]:
    p = line.split('\t')
    V[p[0]] = p[1]
S_pos = int(V['S_pos_0based']); B = 38719; ATTLEN = 146; PL2 = 77437
j_att_mid = S_pos + ATTLEN // 2


def run(cmd, out=None):
    LOG.write('CMD: %s\n' % ' '.join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True, encoding='utf-8', errors='replace')
    LOG.write('  rc=%d stderr: %s\n' % (r.returncode, (r.stderr or '').strip()[:200]))
    if out is not None:
        with open(out, 'w', encoding='utf-8') as f:
            f.write(r.stdout or '')
    if r.returncode != 0:
        raise SystemExit('samtools failed: %s' % ' '.join(cmd))
    return r.stdout


for tag in ['full', 'panel']:
    sam = os.path.join(MAP, '%s_db.sam.gz' % tag)
    if not os.path.exists(sam):
        print('missing', sam); continue
    tmp = os.path.join(MAP, '%s_db.tmp.bam' % tag)
    srt = os.path.join(MAP, '%s_db.sorted.bam' % tag)
    run([ST, 'view', '-b', '-o', tmp, sam])
    run([ST, 'sort', '-m', '1G', '-@', '4', '-o', srt, tmp])
    run([ST, 'index', srt])
    run([ST, 'flagstat', srt], out=os.path.join(MAP, '%s_db.flagstat.txt' % tag))
    run([ST, 'idxstats', srt], out=os.path.join(MAP, '%s_db.idxstats.txt' % tag))
    os.remove(tmp)
    print('done', tag, flush=True)

# junction ±2kb BAM 子集
sub = os.path.join(MAP, 'bam_subsets')
os.makedirs(sub, exist_ok=True)
panel_bam = os.path.join(MAP, 'panel_db.sorted.bam')
regs = [('attJ_junction:3001-7000'), ('retained_window:3074-7073'),
        ('unit1_to_unit2_junction:3001-7000'), ('unit2_to_unit1_circular_junction:3001-7000'),
        ('neg_attJ_shuffled:3001-7000'), ('neg_ctrl_chr08_1:3001-7000')]
for reg in regs:
    name = reg.split(':')[0]
    run([ST, 'view', '-b', '-o', os.path.join(sub, '%s_center4kb.bam' % name), panel_bam, reg])

# samtools depth（算法A）关键区域
full_bam = os.path.join(MAP, 'full_db.sorted.bam')
regions = [
    ('chr08_whole', ['chr08_full:1-2638187']),
    ('pl2_whole', ['pl2_asassembled:1-77437']),
    ('pl1_whole', ['pl1_smbu08:1-128837']),
    ('attJ_2kb', ['chr08_full:%d-%d' % (j_att_mid - 2000 + 1, j_att_mid + 2000)]),
    ('retained_X_internal', ['chr06_full:%d-%d' % (S_pos + ATTLEN + 1, S_pos + B)]),
    ('dimer_forward_2kb', ['pl2_asassembled:%d-%d' % (B - 2000 + 1, B + 2000)]),
    ('dimer_back_wrap_2kb', ['pl2_asassembled:%d-%d' % (PL2 - 2000 + 1, PL2), 'pl2_asassembled:1-2000']),
    ('monomer_mid_2kb', ['pl2_asassembled:%d-%d' % (B // 2 - 2000 + 1, B // 2 + 2000)]),
    ('attL_ortholog_2kb_chr06', ['chr06_full:%d-%d' % (S_pos - 2000 + 1, S_pos + 2000)]),
]
rows = []
for name, regs_ in regions:
    vals = []
    for reg in regs_:
        r = subprocess.run([ST, 'depth', '-a', '-Q', '20', '-r', reg, full_bam],
                           capture_output=True, text=True, encoding='utf-8', errors='replace')
        for line in r.stdout.splitlines():
            c = line.split('\t')
            if len(c) >= 3:
                vals.append(int(c[2]))
    if vals:
        vals_sorted = sorted(vals)
        n = len(vals)
        mean = sum(vals) / n
        med = vals_sorted[n // 2]
        q1 = vals_sorted[n // 4]; q3 = vals_sorted[3 * n // 4]
        rows.append([name, ';'.join(regs_), n, round(mean, 1), med, '%d-%d' % (q1, q3)])
        print(name, rows[-1], flush=True)
with open(os.path.join(OUT, 'depth_samtools_regions.tsv'), 'w', encoding='utf-8') as f:
    f.write('region\tregions(samtools)\tn_positions\tdepth_mean\tdepth_median\tdepth_IQR\n')
    for r_ in rows:
        f.write('\t'.join(str(x) for x in r_) + '\n')

print('bam_convert_depth done')
LOG.write('DONE\n')
