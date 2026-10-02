# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
build_docx.py — Markdown → DOCX（标题/段落/列表/表格；Arial+宋体回退）
输入：manuscript_full_v1.md, 中文完整分析报告.md, figure_legends.md
输出：manuscript_full_v1.docx, 中文完整分析报告.docx, figure_legends.docx
"""
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
from p3lib import WORK
from docx import Document
from docx.shared import Pt, RGBColor
from docx.oxml.ns import qn

MAN = os.path.join(WORK, '11_manuscript')

BOLD_RE = re.compile(r'\*\*(.+?)\*\*')


def add_runs(par, text):
    """处理 **bold** 与 *italic* 简单内联"""
    pos = 0
    for m in BOLD_RE.finditer(text):
        if m.start() > pos:
            par.add_run(text[pos:m.start()])
        r = par.add_run(m.group(1))
        r.bold = True
        pos = m.end()
    if pos < len(text):
        par.add_run(text[pos:])


def md_to_docx(md_path, docx_path, title=None):
    doc = Document()
    style = doc.styles['Normal']
    style.font.name = 'Arial'
    style.font.size = Pt(10.5)
    style._element.rPr.rFonts.set(qn('w:eastAsia'), '宋体')

    lines = open(md_path, encoding='utf-8').read().split('\n')
    i = 0
    while i < len(lines):
        line = lines[i].rstrip()
        if not line.strip():
            i += 1
            continue
        # 表格
        if line.startswith('|') and i + 1 < len(lines) and re.match(r'^\|[\s\-:|]+\|$', lines[i + 1].strip()):
            header = [c.strip() for c in line.strip('|').split('|')]
            rows = []
            j = i + 2
            while j < len(lines) and lines[j].strip().startswith('|'):
                rows.append([c.strip() for c in lines[j].strip('|').split('|')])
                j += 1
            t = doc.add_table(rows=1 + len(rows), cols=len(header))
            t.style = 'Light Grid Accent 1'
            for c, h in enumerate(header):
                cell = t.cell(0, c)
                cell.text = ''
                add_runs(cell.paragraphs[0], h)
            for r, row in enumerate(rows, start=1):
                for c in range(len(header)):
                    cell = t.cell(r, c)
                    cell.text = ''
                    add_runs(cell.paragraphs[0], row[c] if c < len(row) else '')
            for row in t.rows:
                for cell in row.cells:
                    for par in cell.paragraphs:
                        for run in par.runs:
                            run.font.size = Pt(8)
            i = j
            continue
        if line.startswith('# '):
            h = doc.add_heading('', level=0)
            add_runs(h, line[2:])
        elif line.startswith('## '):
            h = doc.add_heading('', level=1)
            add_runs(h, line[3:])
        elif line.startswith('### '):
            h = doc.add_heading('', level=2)
            add_runs(h, line[4:])
        elif line.startswith('- ') or line.startswith('* '):
            par = doc.add_paragraph(style='List Bullet')
            add_runs(par, line[2:])
        elif re.match(r'^\d+\. ', line):
            par = doc.add_paragraph(style='List Number')
            add_runs(par, re.sub(r'^\d+\. ', '', line))
        else:
            par = doc.add_paragraph()
            add_runs(par, line)
        i += 1
    doc.save(docx_path)
    print('saved', docx_path)


md_to_docx(os.path.join(MAN, 'manuscript_full_v1.md'), os.path.join(MAN, 'manuscript_full_v1.docx'))
md_to_docx(os.path.join(MAN, '中文完整分析报告.md'), os.path.join(MAN, '中文完整分析报告.docx'))
md_to_docx(os.path.join(MAN, 'figure_legends.md'), os.path.join(MAN, 'figure_legends.docx'))
print('build_docx done')
