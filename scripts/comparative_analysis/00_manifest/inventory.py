# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
inventory.py — 输入文件清单与 SHA-256（模块 00）
枚举：
  A. smbu06 交付目录全部文件
  B. smbu08 交付目录全部文件
  C. 关键参考文件（论文初稿、创新框架、参考模板、前期工作区关键表）
输出：00_manifest/input_inventory.tsv
      列：sample / category / rel_path / abs_path / size_bytes / mtime / sha256
用法：python inventory.py [--fast]   --fast 跳过 >200MB 文件的哈希（仍记录大小）
"""
import os
import sys
import time
import datetime

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p3lib import WORK, IN06, IN08, P3WS, sha256_file

OUT = os.path.join(WORK, '00_manifest', 'input_inventory.tsv')

CATEGORIES = {
    '1.Cleandata': 'reads',
    '2.Assembly': 'assembly',
    '3.Genome_Component': 'components',
    '4.Genome_Function': 'function',
    '6.Circles_Graphs': 'circos',
}

TARGETS = []


def add_tree(sample, root):
    for dirpath, dirnames, filenames in os.walk(root):
        for fn in filenames:
            p = os.path.join(dirpath, fn)
            rel = os.path.relpath(p, root)
            top = rel.split(os.sep)[0]
            cat = CATEGORIES.get(top, 'root')
            TARGETS.append((sample, cat, rel, p, os.path.getsize(p), os.path.getmtime(p)))


def add_single(sample, cat, p, label=None):
    if os.path.exists(p):
        TARGETS.append((sample, cat, label or os.path.basename(p), p,
                        os.path.getsize(p), os.path.getmtime(p)))
    else:
        TARGETS.append((sample, cat, (label or os.path.basename(p)), p, -1, 0))


def main():
    fast = '--fast' in sys.argv
    add_tree('smbu06', IN06)
    add_tree('smbu08', IN08)
    # 关键参考
    add_single('reference', 'draft', P3WS + r'\smbu06_smbu08_比较基因组论文初稿_v0.1.docx')
    add_single('reference', 'draft', P3WS + r'\smbu06_smbu08_比较研究创新框架.md')
    add_single('reference', 'draft', P3WS + r'\新建 DOCX 文档.docx')
    add_single('reference', 'task', P3WS + r'\task_prompt_smbu06_smbu08_全自动计算分析与正文任务.md')
    add_single('reference', 'task', P3WS + r'\task_prompt_F119_分析任务提示词.md')
    # 参考目录（模板与格式参考）——只取顶层文件
    refdir = P3WS + r'\参考'
    if os.path.isdir(refdir):
        for fn in sorted(os.listdir(refdir)):
            p = os.path.join(refdir, fn)
            if os.path.isfile(p):
                add_single('reference', 'template', p)
    # 前期工作区关键报告
    prev = [
        (P3WS + r'\分析工作区\09_report\中文分析报告.md', 'prev_report'),
        (P3WS + r'\分析工作区\09_report\English_Results_Discussion_draft.md', 'prev_report'),
        (P3WS + r'\分析工作区\09_report\方法与阈值说明.md', 'prev_report'),
        (P3WS + r'\F119_分析工作区\01_audit\菌株身份核验与提交前技术核查_2026-10-01.md', 'prev_report'),
        (P3WS + r'\F119_分析工作区\01_audit\ani_results.tsv', 'prev_table'),
    ]
    for p, cat in prev:
        add_single('previous_work', cat, p)
    # 面板元数据（前期收集的accession列表，用于新面板重建）
    add_single('previous_work', 'prev_table', P3WS + r'\分析工作区\02_panel\accessions.txt')
    add_single('previous_work', 'prev_table', P3WS + r'\分析工作区\02_panel\panel_metadata.tsv')

    t0 = time.time()
    n_hashed = 0
    with open(OUT, 'w', encoding='utf-8') as f:
        f.write('sample\tcategory\trel_path\tabs_path\tsize_bytes\tmtime\tsha256\n')
        for sample, cat, rel, p, size, mt in TARGETS:
            mtime = datetime.datetime.fromtimestamp(mt).strftime('%Y-%m-%d %H:%M:%S') if mt else ''
            do_hash = os.path.exists(p)
            if fast and size > 200 * 1024 * 1024:
                do_hash = False
            h = ''
            if do_hash:
                h = sha256_file(p)
                n_hashed += 1
                if n_hashed % 100 == 0:
                    print('  hashed %d files... (%.0fs)' % (n_hashed, time.time() - t0), flush=True)
            elif size > 0:
                h = 'SKIPPED_FAST_MODE' if fast else ''
            f.write('\t'.join([sample, cat, rel, p, str(size), mtime, h]) + '\n')
    print('wrote %s  (%d entries, %d hashed, %.0fs)' % (OUT, len(TARGETS), n_hashed, time.time() - t0))


if __name__ == '__main__':
    main()
