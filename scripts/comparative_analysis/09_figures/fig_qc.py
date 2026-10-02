# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
fig_qc.py — 图件程序化视觉检查（PyMuPDF）
检查：文本是否越界（超出页面内容区）、文本块两两重叠、字体嵌入、页面尺寸、元素计数
用法：python fig_qc.py <pdf> [<pdf> ...] → 输出 QC 报告行
"""
import sys
import fitz

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

MARGIN_PT = 4.0   # 允许的出血量


def overlap(a, b):
    x0 = max(a[0], b[0]); y0 = max(a[1], b[1])
    x1 = min(a[2], b[2]); y1 = min(a[3], b[3])
    if x1 <= x0 or y1 <= y0:
        return 0.0
    inter = (x1 - x0) * (y1 - y0)
    amin = min((a[2] - a[0]) * (a[3] - a[1]), (b[2] - b[0]) * (b[3] - b[1]))
    return inter / amin if amin else 0.0


def qc(pdf):
    doc = fitz.open(pdf)
    issues = []
    for pno, page in enumerate(doc):
        W, H = page.rect.width, page.rect.height
        d = page.get_text('dict')
        spans = []
        for block in d.get('blocks', []):
            if block.get('type') != 0:
                continue
            for line in block.get('lines', []):
                for sp in line.get('spans', []):
                    if sp['text'].strip():
                        spans.append((sp['bbox'], sp['text']))
        for (x0, y0, x1, y1), t in spans:
            if x0 < -MARGIN_PT or y0 < -MARGIN_PT or x1 > W + MARGIN_PT or y1 > H + MARGIN_PT:
                issues.append('P%d text out of page: %r @(%.1f,%.1f,%.1f,%.1f) page=%.0fx%.0f'
                              % (pno + 1, t[:40], x0, y0, x1, y1, W, H))
        for i in range(len(spans)):
            for j in range(i + 1, len(spans)):
                ov = overlap(spans[i][0], spans[j][0])
                if ov > 0.55:
                    # 同文本（重复绘制）不报
                    if spans[i][1].strip() == spans[j][1].strip():
                        continue
                    issues.append('P%d overlap %.0f%%: %r <> %r'
                                  % (pno + 1, ov * 100, spans[i][1][:38], spans[j][1][:38]))
        fonts = page.get_fonts()
        non_embedded = [f for f in fonts if not f[3]]
        if non_embedded:
            issues.append('P%d non-embedded fonts: %s' % (pno + 1, [f[3] for f in non_embedded]))
    print('=== QC %s | pages=%d | spans=%d | issues=%d' % (pdf, len(doc), len(spans), len(issues)))
    for it in issues[:30]:
        print('   ', it)
    return len(issues)


if __name__ == '__main__':
    total = 0
    for p in sys.argv[1:]:
        total += qc(p)
    print('TOTAL ISSUES:', total)
