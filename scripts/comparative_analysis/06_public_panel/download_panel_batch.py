# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
download_panel_batch.py — 模块 06：面板批量下载（替代逐条模式；每批 25 株单次 datasets 调用）
断点续传：已存在 .fna 的 accession 跳过；批 zip 下载→解包到 <acc>/→删 zip
"""
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
from p3lib import WORK, DATASETS, open_log

OUT = os.path.join(WORK, '06_public_panel')
DL = os.path.join(OUT, 'downloads')
os.makedirs(DL, exist_ok=True)
LOG = open_log(os.path.join(OUT, 'logs', 'download_panel.log'))

with open(os.path.join(OUT, 'metadata', 'panel_accessions.txt'), encoding='utf-8') as f:
    accs = [l.strip() for l in f if l.strip()]


def has_fna(acc):
    d = os.path.join(DL, acc)
    if not os.path.isdir(d):
        return None
    for root, _, files in os.walk(d):
        for fn in files:
            if fn.endswith('.fna'):
                return os.path.join(root, fn)
    return None


todo = [a for a in accs if not has_fna(a)]
LOG.write('total=%d todo=%d\n' % (len(accs), len(todo)))
BATCH = 25
for bi in range(0, len(todo), BATCH):
    batch = todo[bi:bi + BATCH]
    zipf = os.path.join(DL, '_batch_%d.zip' % (bi // BATCH))
    cmd = [DATASETS, 'download', 'genome', 'accession'] + batch + \
          ['--include', 'genome', '--filename', zipf]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    if r.returncode != 0 or not os.path.exists(zipf):
        LOG.write('batch %d FAILED: %s\n' % (bi // BATCH, (r.stderr or '')[-300:]))
        continue
    # 解包到临时目录后按 accession 归位
    tmpd = os.path.join(DL, '_batch_%d' % (bi // BATCH))
    os.makedirs(tmpd, exist_ok=True)
    subprocess.run(['unzip', '-o', '-q', zipf, '-d', tmpd], timeout=1800)
    # 结构: ncbi_dataset/data/<acc>/<acc>_*.fna
    data = os.path.join(tmpd, 'ncbi_dataset', 'data')
    n_ok = 0
    if os.path.isdir(data):
        # 处理合并模式（多 accession 时目录按各自的 acc）
        for entry in os.listdir(data):
            src = os.path.join(data, entry)
            if not os.path.isdir(src):
                continue
            if entry in batch:
                dst = os.path.join(DL, entry)
                os.makedirs(dst, exist_ok=True)
                shutil.copytree(src, os.path.join(dst, 'data'), dirs_exist_ok=True)
                n_ok += 1
        # 多 accession 下载可能是扁平结构，用文件名匹配
        for entry in os.listdir(data):
            if os.path.isdir(os.path.join(data, entry)):
                continue
            for acc in batch:
                if entry.startswith(acc):
                    dst = os.path.join(DL, acc, 'data', entry)
                    os.makedirs(os.path.dirname(dst), exist_ok=True)
                    shutil.copy2(os.path.join(data, entry), dst)
                    n_ok += 1
    LOG.write('batch %d: %s -> %d placed\n' % (bi // BATCH, batch, n_ok))
    try:
        os.remove(zipf)
        shutil.rmtree(tmpd)
    except OSError:
        pass
    print('batch %d done (%d)' % (bi // BATCH, n_ok), flush=True)

# 汇总
rows = []
for a in accs:
    fna = has_fna(a)
    rows.append((a, 'ok' if fna else 'MISSING', fna or ''))
with open(os.path.join(OUT, 'metadata', 'download_log.tsv'), 'w', encoding='utf-8') as f:
    f.write('accession\tstatus\tpath\n')
    for r in rows:
        f.write('\t'.join(r) + '\n')
n_ok = sum(1 for r in rows if r[1] == 'ok')
print('download_panel_batch done: %d/%d' % (n_ok, len(rows)))
LOG.write('DONE: %d/%d ok\n' % (n_ok, len(rows)))
