# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""retry_missing.py — 重试下载失败的面板条目（单条模式，记录错误原因）"""
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
from p3lib import WORK, DATASETS, open_log

OUT = os.path.join(WORK, '06_public_panel')
DL = os.path.join(OUT, 'downloads')
LOG = open_log(os.path.join(OUT, 'logs', 'download_panel.log'))

missing = []
with open(os.path.join(OUT, 'metadata', 'download_log.tsv'), encoding='utf-8') as f:
    next(f)
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) >= 2 and p[1] == 'MISSING':
            missing.append(p[0])
LOG.write('retry missing: %d\n' % len(missing))
for a in missing:
    zipf = os.path.join(DL, a + '.zip')
    r = subprocess.run([DATASETS, 'download', 'genome', 'accession', a, '--include', 'genome',
                        '--filename', zipf], capture_output=True, text=True, timeout=600)
    if r.returncode != 0 or not os.path.exists(zipf):
        LOG.write('RETRY-FAIL %s: %s\n' % (a, (r.stderr or r.stdout or '')[-200:].replace('\n', ' ')))
        continue
    d = os.path.join(DL, a)
    os.makedirs(d, exist_ok=True)
    subprocess.run(['unzip', '-o', '-q', zipf, '-d', d], timeout=600)
    try:
        os.remove(zipf)
    except OSError:
        pass
    LOG.write('RETRY-OK %s\n' % a)
    print('ok', a, flush=True)

# 重新汇总
accs = [l.strip() for l in open(os.path.join(OUT, 'metadata', 'panel_accessions.txt'), encoding='utf-8') if l.strip()]
rows = []
for a in accs:
    fna = None
    base = os.path.join(DL, a)
    if os.path.isdir(base):
        for root, _, files in os.walk(base):
            for fn in files:
                if fn.endswith('.fna'):
                    fna = os.path.join(root, fn)
    rows.append((a, 'ok' if fna else 'MISSING', fna or ''))
with open(os.path.join(OUT, 'metadata', 'download_log.tsv'), 'w', encoding='utf-8') as f:
    f.write('accession\tstatus\tpath\n')
    for r in rows:
        f.write('\t'.join(r) + '\n')
n_ok = sum(1 for r in rows if r[1] == 'ok')
print('after retry: %d/%d ok' % (n_ok, len(rows)))
LOG.write('after retry: %d/%d ok\n' % (n_ok, len(rows)))
