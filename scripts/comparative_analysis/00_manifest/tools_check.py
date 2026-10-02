# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""tools_check.py — 记录软件与系统版本（模块 00）。输出 00_manifest/tool_versions.txt"""
import os
import platform
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p3lib import (WORK, BLAST, MM2, MAFFT, IQTREE, PY, DATASETS)

OUT = os.path.join(WORK, '00_manifest', 'tool_versions.txt')

CMDS = [
    ('python', [PY, '--version']),
    ('minimap2', [MM2, '--version']),
    ('blastn', [BLAST + r'\blastn.exe', '-version']),
    ('makeblastdb', [BLAST + r'\makeblastdb.exe', '-version']),
    ('tblastn', [BLAST + r'\tblastn.exe', '-version']),
    ('blastp', [BLAST + r'\blastp.exe', '-version']),
    ('mafft', ['cmd', '/c', MAFFT, '--version']),
    ('iqtree2', [IQTREE, '--version']),
    ('datasets', [DATASETS, 'version']),
    ('java', ['java', '-version']),
    ('Rscript', ['Rscript', '--version']),
    ('xelatex', ['xelatex', '--version']),
]

PYMODS = ['numpy', 'pandas', 'scipy', 'matplotlib', 'docx', 'openpyxl', 'pyrodigal', 'multiqc', 'PIL']


def main():
    lines = []
    lines.append('smbu06_smbu08 比较基因组完整分析 — 软件版本记录')
    lines.append('生成时间: %s' % subprocess.run([PY, '-c', 'import datetime;print(datetime.datetime.now())'],
                                                 capture_output=True, text=True).stdout.strip())
    lines.append('平台: %s' % platform.platform())
    lines.append('---')
    for name, cmd in CMDS:
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=60,
                               encoding='utf-8', errors='replace')
            out = (r.stdout or '').strip() or (r.stderr or '').strip()
            out = out.splitlines()[0] if out else ''
            lines.append('%s: %s' % (name, out))
        except Exception as e:
            lines.append('%s: ERROR %s' % (name, e))
    lines.append('---python modules---')
    code = ';'.join(
        "import importlib;m=importlib.import_module('%s');print('%s',getattr(m,'__version__','?'))" % (m, m)
        for m in PYMODS)
    r = subprocess.run([PY, '-c', code], capture_output=True, text=True, encoding='utf-8', errors='replace')
    lines.append(r.stdout.strip())
    lines.append('---说明---')
    lines.append('Windows 环境无 WSL/Docker/编译器；minimap2 使用社区 MSYS2 构建 (win-ngs/minimap2-windows-build v2.31-r1302)；')
    lines.append('BLAST+/MAFFT/IQ-TREE/datasets 为官方 Windows 构建；中文路径经 ASCII junction 中转')
    lines.append('(<PROJECT_ROOT>\\{smbu06,smbu08,p3cmp})，因 BLAST+ LMDB 不支持非 ASCII 路径。')
    with open(OUT, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')
    print('\n'.join(lines))
    print('wrote', OUT)


if __name__ == '__main__':
    main()
