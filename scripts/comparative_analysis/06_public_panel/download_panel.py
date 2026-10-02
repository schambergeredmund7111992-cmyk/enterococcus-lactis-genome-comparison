# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
download_panel.py — 模块 06：公共基因组面板下载（可断点续传）
来源：前期面板 accession 清单（smbu06 分析工作区，208 株，来自 NCBI Datasets）
      + nisin 阳性功能对照（前期已验证 5 株）+ 本研究参考株（194/14B4/MG1363/IL1403/F44 区段）
输出：06_public_panel/downloads/<accession>/ 及 download_log.tsv
"""
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
from p3lib import WORK, DATASETS, open_log

OUT = os.path.join(WORK, '06_public_panel')
DL = os.path.join(OUT, 'downloads')
os.makedirs(DL, exist_ok=True)
LOG = open_log(os.path.join(OUT, 'logs', 'download_panel.log'))

accs = []
seen = set()
src = WORK + r'\..\p3ws\分析工作区\02_panel\accessions.txt'
src = r'<PROJECT_ROOT>\分析工作区\02_panel\accessions.txt'
with open(src, encoding='utf-8') as f:
    for line in f:
        a = line.strip()
        if a and a not in seen:
            seen.add(a)
            accs.append(a)
LOG.write('base accessions: %d\n' % len(accs))

NISIN_POS = ['GCA_000344575.1', 'GCA_000761115.1', 'GCA_000807375.1', 'GCA_002804185.1', 'GCA_002804285.1']
for a in NISIN_POS:
    if a not in seen:
        seen.add(a)
        accs.append(a)
LOG.write('with nisin-positive controls: %d\n' % len(accs))

with open(os.path.join(OUT, 'metadata', 'panel_accessions.txt'), 'w', encoding='utf-8') as f:
    f.write('\n'.join(accs) + '\n')

rows = []
for i, a in enumerate(accs):
    d = os.path.join(DL, a)
    fna = None
    if os.path.isdir(d):
        for root, _, files in os.walk(d):
            for fn in files:
                if fn.endswith('.fna'):
                    fna = os.path.join(root, fn)
                    break
            if fna:
                break
    if fna:
        rows.append((a, 'exists', fna))
        continue
    zipf = os.path.join(DL, a + '.zip')
    cmd = [DATASETS, 'download', 'genome', 'accession', a, '--include', 'genome', '--filename', zipf]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    ok = r.returncode == 0
    if ok:
        os.makedirs(d, exist_ok=True)
        subprocess.run(['unzip', '-o', '-q', zipf, '-d', d], timeout=300)
        try:
            os.remove(zipf)
        except OSError:
            pass
        for root, _, files in os.walk(d):
            for fn in files:
                if fn.endswith('.fna'):
                    fna = os.path.join(root, fn)
        rows.append((a, 'downloaded', fna or ''))
    else:
        rows.append((a, 'FAILED', (r.stderr or '')[-200:].replace('\n', ' ')))
    LOG.write('[%d/%d] %s: %s\n' % (i + 1, len(accs), a, rows[-1][1]))
    if (i + 1) % 20 == 0:
        print('progress %d/%d' % (i + 1, len(accs)), flush=True)

with open(os.path.join(OUT, 'metadata', 'download_log.tsv'), 'w', encoding='utf-8') as f:
    f.write('accession\tstatus\tpath_or_error\n')
    for r in rows:
        f.write('\t'.join(r) + '\n')
n_ok = sum(1 for r in rows if r[1] in ('exists', 'downloaded'))
LOG.write('DONE: %d/%d ok\n' % (n_ok, len(rows)))
print('download_panel done: %d/%d ok' % (n_ok, len(rows)))
