# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
build_docx_v2.py — Markdown → DOCX（含图片内嵌、图注、上标、表格）
支持语法：
  [[FIG:相对路径|宽cm]]    ← 插入图片（居中，附于上一段图注之后）
  ^text^                   ← 上标（作者/机构标记）
  **bold**                 ← 粗体
  | 表格 |                 ← 表格
输出：manuscript_full_v2.docx（英文，图内嵌）、中文完整分析报告_v2.docx、figure_legends.docx
"""
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
from p3lib import WORK
from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn

MAN = os.path.join(WORK, '11_manuscript')
BOLD_RE = re.compile(r'(\*\*.+?\*\*)')
SUP_RE = re.compile(r'(\^[^\^]+\^)')
FIG_RE = re.compile(r'^\[\[FIG:([^|\]]+)\|([\d.]+)\]\]$')


def add_runs(par, text):
    for chunk in BOLD_RE.split(text):
        if not chunk:
            continue
        if chunk.startswith('**') and chunk.endswith('**'):
            r = par.add_run(chunk[2:-2])
            r.bold = True
            continue
        for piece in SUP_RE.split(chunk):
            if not piece:
                continue
            if piece.startswith('^') and piece.endswith('^') and len(piece) > 2:
                r = par.add_run(piece[1:-1])
                r.font.superscript = True
            else:
                par.add_run(piece)


def md_to_docx(md_path, docx_path):
    doc = Document()
    style = doc.styles['Normal']
    style.font.name = 'Arial'
    style.font.size = Pt(10.5)
    style._element.rPr.rFonts.set(qn('w:eastAsia'), '宋体')
    lines = open(md_path, encoding='utf-8').read().split('\n')
    n_img = 0
    i = 0
    while i < len(lines):
        line = lines[i].rstrip()
        if not line.strip():
            i += 1
            continue
        mfig = FIG_RE.match(line.strip())
        if mfig:
            rel, width = mfig.group(1), float(mfig.group(2))
            img = os.path.join(WORK, rel.replace('/', os.sep))
            if os.path.exists(img):
                p = doc.add_paragraph()
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER
                p.add_run().add_picture(img, width=Cm(width))
                n_img += 1
            else:
                doc.add_paragraph('[missing figure: %s]' % rel)
            i += 1
            continue
        if line.strip() == '---':
            i += 1
            continue
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
            for r_, row in enumerate(rows, start=1):
                for c in range(len(header)):
                    cell = t.cell(r_, c)
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
        elif line.startswith('Figure. ') or line.startswith('图 '):
            par = doc.add_paragraph()
            r = par.add_run(line)
            r.font.size = Pt(9)
            r.italic = True
        else:
            par = doc.add_paragraph()
            add_runs(par, line)
        i += 1
    doc.save(docx_path)
    print('saved %s (images embedded: %d)' % (docx_path, n_img))


md_to_docx(os.path.join(MAN, 'manuscript_full_v2.md'), os.path.join(MAN, 'manuscript_full_v2.docx'))
if os.path.exists(os.path.join(MAN, '中文完整分析报告_v2.md')):
    md_to_docx(os.path.join(MAN, '中文完整分析报告_v2.md'), os.path.join(MAN, '中文完整分析报告_v2.docx'))
if os.path.exists(os.path.join(MAN, 'figure_legends_v2.md')):
    md_to_docx(os.path.join(MAN, 'figure_legends_v2.md'), os.path.join(MAN, 'figure_legends_v2.docx'))
print('build_docx_v2 done')
