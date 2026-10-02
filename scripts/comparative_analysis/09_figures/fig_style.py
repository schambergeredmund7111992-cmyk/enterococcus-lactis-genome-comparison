# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""fig_style.py — 统一图件风格（Nature 风格单栏/双栏；英文标注）"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Rectangle, FancyArrowPatch

plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 7,
    'axes.linewidth': 0.6,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'xtick.major.size': 2.5,
    'ytick.major.size': 2.5,
    'legend.frameon': False,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'svg.fonttype': 'none',
})

# 统一配色
C = {
    'smbu06': '#2166AC',
    'smbu08': '#B2182B',
    'shared': '#4D4D4D',
    'prophage': '#E08214',
    'plasmid2': '#8073AC',
    'bacteriocin': '#D53E4F',
    'immunity': '#F46D43',
    'mge': '#999999',
    'core': '#66C2A5',
    'other': '#BFBFBF',
    'highlight': '#FDAE61',
}


def savefig(fig, base):
    """导出 PNG(600dpi)/PDF/SVG 三格式；base 为不含扩展名的路径"""
    import os
    os.makedirs(os.path.dirname(base), exist_ok=True)
    fig.savefig(base + '.png', dpi=600, bbox_inches='tight')
    fig.savefig(base + '.pdf', bbox_inches='tight')
    fig.savefig(base + '.svg', bbox_inches='tight')
    print('saved', base + '.{png,pdf,svg}')


def gene_arrow(ax, x0, x1, y, strand, color, height=0.32, label=None, fontsize=5.2,
               lw=0.4, label_below=False, text_color='black'):
    """基因箭头（坐标以 kb 计）"""
    x0, x1 = min(x0, x1), max(x0, x1)
    w = max(x1 - x0, 0.15)
    # 箭头主体
    if strand == '+':
        body = Rectangle((x0, y - height / 2), w - min(w * 0.35, 0.9), height,
                         facecolor=color, edgecolor='black', linewidth=lw, zorder=3)
        tip = FancyArrow(x0 + w - min(w * 0.35, 0.9), y, min(w * 0.35, 0.9), 0,
                         width=height, head_width=height * 1.35, head_length=min(w * 0.35, 0.9),
                         length_includes_head=True, facecolor=color, edgecolor='black',
                         linewidth=lw, zorder=3)
    else:
        body = Rectangle((x0 + min(w * 0.35, 0.9), y - height / 2), w - min(w * 0.35, 0.9), height,
                         facecolor=color, edgecolor='black', linewidth=lw, zorder=3)
        tip = FancyArrow(x0 + min(w * 0.35, 0.9), y, -min(w * 0.35, 0.9), 0,
                         width=height, head_width=height * 1.35, head_length=min(w * 0.35, 0.9),
                         length_includes_head=True, facecolor=color, edgecolor='black',
                         linewidth=lw, zorder=3)
    ax.add_patch(body)
    ax.add_patch(tip)
    if label:
        yy = y - height * 0.95 if label_below else y + height * 0.95
        ax.text((x0 + x1) / 2, yy, label, ha='center', va='top' if label_below else 'bottom',
                fontsize=fontsize, color=text_color, zorder=4)


class PanelLabel:
    def __init__(self, fig, x=-0.02, y=1.02):
        self.fig = fig
        self.x, self.y = x, y

    def __call__(self, ax, letter, dx=0.0, dy=0.0, fontsize=9):
        ax.text(self.x + dx, self.y + dy, letter, transform=ax.transAxes, fontsize=fontsize,
                fontweight='bold', va='bottom', ha='right')
