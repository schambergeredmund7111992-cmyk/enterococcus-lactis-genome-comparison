# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
build_docx_final.py — 严格模仿 参考\Manuscript.docx 排版生成最终 DOCX
排版规范（从参考文档 programmatically 提取）：
  页面 A4；边距 L/R 3.17cm，T/B 2.54cm
  正文 Times New Roman 10.5pt（w:sz=21），行距 2.0（w:line=480 auto），首行缩进 18pt（360 twips），段前后 0
  标题 1/2/3：Times New Roman 12pt 粗体；1 级/2 级为 “N.<TAB>题名”，3 级为 “N.N.N 题名”
  Title：TNR 16pt 粗体居中；References: 居中 TNR 20pt；参考文献条目 TNR 10pt
  图注 “Figure. N …” 整行粗体、无缩进、双倍行距；图片居中、置于题注上方（与参考文档一致）
  表题 “Table.N …” 粗体居中
输出：manuscript_final.docx（图内嵌）、中文完整分析报告_final.docx
"""
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
from p3lib import WORK
from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

MAN = os.path.join(WORK, '11_manuscript')
BOLD_RE = re.compile(r'(\*\*.+?\*\*)')
SUP_RE = re.compile(r'(\^[^\^]+\^)')
FIG_RE = re.compile(r'^\[\[FIG:([^|\]]+)\|([\d.]+)\]\]$')
H1_RE = re.compile(r'^(\d+\.[\s]?[A-Z].*|[A-Z][a-z].*)$')   # 由专用列表控制
H2_RE = re.compile(r'^\d+\.\d+[\s]')
H3_RE = re.compile(r'^\d+\.\d+\.\d+[\s]')


def set_font(run, size=None, bold=None, italic=None, east=None):
    run.font.name = 'Times New Roman'
    rPr = run._element.get_or_add_rPr()
    rf = rPr.find(qn('w:rFonts'))
    if rf is None:
        rf = OxmlElement('w:rFonts')
        rPr.insert(0, rf)
    rf.set(qn('w:ascii'), 'Times New Roman')
    rf.set(qn('w:hAnsi'), 'Times New Roman')
    rf.set(qn('w:cs'), 'Times New Roman')
    rf.set(qn('w:eastAsia'), east or '宋体')
    if size is not None:
        run.font.size = Pt(size)
        # 同时写 cs 字号，模拟参考文档
        szCs = OxmlElement('w:szCs')
        szCs.set(qn('w:val'), str(int(size * 2)))
        rPr.append(szCs)
    if bold is not None:
        run.font.bold = bold
    if italic is not None:
        run.font.italic = italic


def add_runs(par, text, size=None, base_bold=None):
    for chunk in BOLD_RE.split(text):
        if not chunk:
            continue
        if chunk.startswith('**') and chunk.endswith('**'):
            r = par.add_run(chunk[2:-2])
            set_font(r, size=size, bold=True)
            continue
        for piece in SUP_RE.split(chunk):
            if not piece:
                continue
            if piece.startswith('^') and piece.endswith('^') and len(piece) > 2:
                r = par.add_run(piece[1:-1])
                set_font(r, size=size, bold=base_bold)
                r.font.superscript = True
            else:
                r = par.add_run(piece)
                set_font(r, size=size, bold=base_bold)


def body_par(doc, indent=True):
    p = doc.add_paragraph()
    pf = p.paragraph_format
    pf.line_spacing = 2.0
    pf.space_before = Pt(0)
    pf.space_after = Pt(0)
    if indent:
        pf.first_line_indent = Pt(18)
    return p


def heading(doc, text, level):
    p = doc.add_paragraph()
    pf = p.paragraph_format
    pf.line_spacing = 2.0
    pf.space_before = Pt(0)
    pf.space_after = Pt(0)
    if level == 1 and not text.startswith(tuple('123456789')):
        # Introduction 等无编号一级标题
        r = p.add_run(text)
        set_font(r, size=12, bold=True)
        return p
    m = re.match(r'^(\d+(?:\.\d+)*)\.?([\s\t]?)(.*)$', text)
    if m and level >= 1:
        num, sep, title = m.group(1), m.group(2), m.group(3)
        if level == 1:
            full = '%s.\t%s' % (num, title)      # 参考格式："2.<TAB>题名"
        elif level == 2:
            full = '%s\t%s' % (num, title)        # 参考格式："2.1<TAB>题名"（无句点）
        else:
            full = '%s %s' % (num, title)         # 参考格式："2.2.1 题名"（空格）
        r = p.add_run(full)
        set_font(r, size=12, bold=True)
        return p
    r = p.add_run(text)
    set_font(r, size=12, bold=True)
    return p


def ref_par(doc, text):
    p = doc.add_paragraph()
    pf = p.paragraph_format
    pf.line_spacing = 1.5
    pf.space_before = Pt(0)
    pf.space_after = Pt(0)
    add_runs(p, text, size=10)
    return p


SEC1 = {'Introduction', 'Supplementary Materials', 'Code availability', 'Funding Declaration',
        'Conflicts of Interest', 'Credit Author Statement', 'Data availability', 'Limitations.'}


def md_to_docx(md_path, docx_path, refs_center=True):
    doc = Document()
    # 页面：A4 + 参考文档边距
    sec = doc.sections[0]
    sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)
    sec.left_margin = sec.right_margin = Cm(3.17)
    sec.top_margin = sec.bottom_margin = Cm(2.54)
    # Normal 样式：TNR 10.5pt
    st = doc.styles['Normal']
    st.font.name = 'Times New Roman'
    st.font.size = Pt(10.5)
    st._element.rPr.rFonts.set(qn('w:eastAsia'), '宋体')

    lines = open(md_path, encoding='utf-8').read().split('\n')
    i = 0
    n_img = 0
    pending_img = None   # 缓存 [[FIG]]：等题注输出后再插图（或先图后题）——参考文档为 图→题注
    while i < len(lines):
        line = lines[i].rstrip()
        i += 1
        if not line.strip() or line.strip() == '---':
            continue
        mfig = FIG_RE.match(line.strip())
        if mfig:
            rel, width = mfig.group(1), float(mfig.group(2))
            img = os.path.join(WORK, rel.replace('/', os.sep))
            # 参考文档布局：图片在题注上方 → 将图片插到上一段（题注）之前
            if os.path.exists(img):
                p = doc.add_paragraph()
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER
                pf = p.paragraph_format
                pf.line_spacing = 2.0
                pf.space_before = Pt(0)
                pf.space_after = Pt(0)
                p.add_run().add_picture(img, width=Cm(width))
                n_img += 1
                # 调整顺序：把图片段落移到最后一段文字（题注）之前
                if len(doc.paragraphs) >= 2:
                    prev = doc.paragraphs[-2]
                    if prev.text.strip().startswith(('Figure.', '图 ')):
                        img_p = doc.paragraphs[-1]._element
                        prev._element.addprevious(img_p)
            i_last = i
            continue
        # 标题识别
        stripped = line.strip()
        if stripped in SEC1 or stripped == 'Introduction':
            heading(doc, stripped, 1)
            continue
        if stripped.startswith('References'):
            p = doc.add_paragraph()
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            r = p.add_run('References:')
            set_font(r, size=20, bold=False)
            continue
        m2 = re.match(r'^(\d+\.\d+)\s+(.*)$', stripped)
        if m2:
            heading(doc, stripped, 2)
            continue
        m3 = re.match(r'^(\d+\.\d+\.\d+)\s+(.*)$', stripped)
        if m3:
            heading(doc, stripped, 3)
            continue
        m1 = re.match(r'^(\d+)\.([A-Z].*)$', stripped)
        if m1:
            heading(doc, stripped, 1)
            continue
        if stripped == 'Abstract:':
            # 下一非空行并入同一段落
            j = i
            while j < len(lines) and not lines[j].strip():
                j += 1
            txt = lines[j].strip() if j < len(lines) else ''
            p = body_par(doc)
            r = p.add_run('Abstract: ')
            set_font(r, size=12, bold=True)
            add_runs(p, txt, size=10.5)
            i = j + 1
            continue
        if stripped.startswith('Keywords:'):
            p = body_par(doc, indent=False)
            r = p.add_run('Keywords: ')
            set_font(r, size=10.5, bold=True, italic=True)
            add_runs(p, stripped[len('Keywords:'):].strip(), size=10.5)
            continue
        if stripped.startswith('Figure. ') or stripped.startswith('Table.') or stripped.startswith('Table '):
            p = body_par(doc, indent=False)
            if stripped.startswith('Table'):
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            add_runs(p, stripped, size=10.5, base_bold=True)
            continue
        if re.match(r'^\d+\. ', stripped) and False:
            continue
        if stripped.startswith('- ') or stripped.startswith('* '):
            p = body_par(doc, indent=False)
            p.paragraph_format.left_indent = Pt(18)
            add_runs(p, stripped[2:], size=10.5)
            continue
        # 作者块（占位识别：含 ^ 上标或 # 或 ✉）
        if stripped.startswith('[') and stripped.endswith(']'):
            p = body_par(doc, indent=False)
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            add_runs(p, stripped, size=10.5)
            continue
        if stripped.startswith('^') or (('#' in stripped or '✉' in stripped) and len(stripped) < 300 and '。' not in stripped and '.' not in stripped[:0]):
            p = body_par(doc, indent=False)
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            add_runs(p, stripped, size=10.5)
            continue
        if re.match(r'^(Fan Ding|#These|✉)', stripped):
            p = body_par(doc, indent=False)
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            add_runs(p, stripped, size=10.5)
            continue
        # 参考文献条目（含 DOI 或以 Author, X. 开头）
        if 'doi.org' in stripped or re.match(r'^[A-Z][a-zA-Z\-]+,\s+[A-Z]\.', stripped):
            ref_par(doc, stripped)
            continue
        # 普通正文
        p = body_par(doc, indent=True)
        add_runs(p, stripped, size=10.5)
    doc.save(docx_path)
    print('saved %s (images=%d)' % (docx_path, n_img))


# 主标题特殊处理：第一行 # 标题 → 16pt 粗体居中
def md_to_docx_full(md_path, docx_path):
    # 预先抽取标题行
    lines = open(md_path, encoding='utf-8').read().split('\n')
    tmp = []
    title_done = False
    for ln in lines:
        if not title_done and ln.startswith('# '):
            continue
        tmp.append(ln)
    tmp_md = os.path.join(MAN, '_tmp_body.md')
    open(tmp_md, 'w', encoding='utf-8').write('\n'.join(tmp))
    # 生成正文
    md_to_docx(tmp_md, docx_path)
    # 在文档开头插入标题
    from docx import Document as D2
    doc = D2(docx_path)
    title = lines[0][2:].strip()
    p = doc.paragraphs[0].insert_paragraph_before()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    pf = p.paragraph_format
    pf.line_spacing = 2.0
    pf.space_before = Pt(12)
    r = p.add_run(title)
    set_font(r, size=16, bold=True)
    doc.save(docx_path)
    print('title inserted:', title[:60])


md_to_docx_full(os.path.join(MAN, 'manuscript_full_v2.md'), os.path.join(MAN, 'manuscript_final.docx'))
if os.path.exists(os.path.join(MAN, '中文完整分析报告_v2.md')):
    md_to_docx_full(os.path.join(MAN, '中文完整分析报告_v2.md'), os.path.join(MAN, '中文完整分析报告_final.docx'))
print('build_docx_final done')
