# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
map_reads2.py — 02_mapping（本会话独立实现）：
  对 (A) full_refs（chr06_full, chr08_full, pl1_smbu08, pl2_asassembled）
     与 (B) panel_refs（13 条 junction/模型/对照参考）
  各运行 minimap2 两次：-c → PAF(.paf.gz)；-a --MD --secondary=yes → SAM(.sam.gz)
  保留原始命令/版本/rc/耗时；未做任何过滤。
"""
import gzip
import hashlib
import os
import subprocess
import sys
import time

ATTJ = r'<PROJECT_ROOT>\longread_validation'
REFDIR = os.path.join(ATTJ, '01_references', 'references')
MAP = os.path.join(ATTJ, '02_mapping')
MM2 = r'<TOOLS_ROOT>\minimap2-win\minimap2-2.31-r1302-windows-x86_64-ucrt64\minimap2.exe'
FQ = r'<PROJECT_ROOT>\smbu08\smbu08\1.Cleandata\smbu08.filtered_reads.fq.gz'
LOG = open(os.path.join(ATTJ, '00_manifest', 'logs', 'map_reads2.log'), 'a', encoding='utf-8')

FULL_REFS = ['chr06_full', 'chr08_full', 'pl1_smbu08', 'pl2_asassembled']
PANEL_REFS = ['attJ_junction', 'retained_window', 'unit1', 'unit2', 'unit1_to_unit2_junction',
              'unit2_to_unit1_circular_junction', 'monomer_circle_repA', 'monomer_circle_repB',
              'dimer_back_internal', 'neg_attJ_shuffled', 'neg_ctrl_chr08_1', 'neg_ctrl_chr08_2',
              'neg_ctrl_chr08_3']


def read_fasta(path):
    recs, name, buf = [], None, []
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    recs.append((name, ''.join(buf)))
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line.strip())
    if name is not None:
        recs.append((name, ''.join(buf)))
    return recs


def build(names, out):
    with open(out, 'w') as fo:
        for n in names:
            seq = read_fasta(os.path.join(REFDIR, n + '.fasta'))[0][1]
            fo.write('>%s\n' % n)
            for i in range(0, len(seq), 60):
                fo.write(seq[i:i + 60] + '\n')
    h = hashlib.md5(open(out, 'rb').read()).hexdigest()
    LOG.write('built %s md5=%s refs=%s\n' % (os.path.basename(out), h, ','.join(names)))
    return h


def sha256_file(p, bs=1 << 22):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        while True:
            b = f.read(bs)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def run(cmd, out_path, gz, tag):
    t0 = time.time()
    LOG.write('CMD %s: %s\n' % (tag, ' '.join(cmd)))
    with open(out_path + ('.tmp' if gz else ''), 'wb') as raw:
        if gz:
            p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            with gzip.open(out_path, 'wb') as gzout:
                while True:
                    chunk = p1.stdout.read(1 << 20)
                    if not chunk:
                        break
                    gzout.write(chunk)
            err = p1.stderr.read().decode('utf-8', 'replace')
            rc = p1.wait()
        else:
            p2 = subprocess.run(cmd, stdout=raw, stderr=subprocess.PIPE)
            rc = p2.returncode
            err = p2.stderr.decode('utf-8', 'replace')
    dt = time.time() - t0
    sz = os.path.getsize(out_path)
    LOG.write('RC %s rc=%d sec=%.0f size=%d stderr_tail=%s\n' % (tag, rc, dt, sz, err.strip()[-300:]))
    print('%s rc=%d %.0fs %.1fMB' % (tag, rc, dt, sz / 1e6), flush=True)
    if rc != 0:
        raise SystemExit('minimap2 failed: %s' % tag)


full_fa = os.path.join(MAP, 'full_refs.fasta')
panel_fa = os.path.join(MAP, 'panel_refs.fasta')
build(FULL_REFS, full_fa)
build(PANEL_REFS, panel_fa)

MM2_COMMON = [MM2, '-x', 'map-ont', '-t', '14']
for tag, ref in [('full', full_fa), ('panel', panel_fa)]:
    run(MM2_COMMON + ['-c', '--secondary=yes', ref, FQ], os.path.join(MAP, '%s_db.paf.gz' % tag), True, '%s_paf' % tag)
    run(MM2_COMMON + ['-a', '--MD', '--secondary=yes', ref, FQ], os.path.join(MAP, '%s_db.sam.gz' % tag), True, '%s_sam' % tag)

LOG.write('DONE map_reads2\n')
print('map_reads2 done')
