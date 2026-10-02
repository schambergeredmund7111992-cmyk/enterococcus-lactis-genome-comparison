# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
build_supplementary.py — 模块 10：补充表工作簿（S1–S6）
输出：10_supplement/tables/Supplementary_Tables.xlsx（每表一个 sheet）+ 各表 tsv 副本
"""
import csv
import os
import shutil
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
from p3lib import WORK
from openpyxl import Workbook
from openpyxl.utils import get_column_letter

W = WORK
OUT = os.path.join(W, '10_supplement', 'tables')
os.makedirs(OUT, exist_ok=True)

TABLES = {
    'TableS1_assembly_stats': [
        W + r'\01_assembly_qc\replicon_stats.tsv',
        W + r'\01_assembly_qc\reads_qc\multiqc_report\multiqc_data\multiqc_fastqc.txt',
        W + r'\01_assembly_qc\mapping\mapping_stats.tsv',
        W + r'\01_assembly_qc\mapping\coverage_deconvolved.tsv',
        W + r'\01_assembly_qc\mapping\coverage_per_replicon.tsv',
        W + r'\01_assembly_qc\mapping\coverage_uniformity.tsv',
        W + r'\01_assembly_qc\mapping\contamination_summary.tsv',
    ],
    'TableS2_candidate_genes': [
        W + r'\05_bacteriocins\search\candidate_table.tsv',
        W + r'\05_bacteriocins\search\physchem_candidates.tsv',
        W + r'\05_bacteriocins\search\neighborhood_GL002641.tsv',
        W + r'\05_bacteriocins\search\neighborhood_GL002679.tsv',
        W + r'\05_bacteriocins\search\alignment_GL002641_vs_enterocinP.txt',
        W + r'\05_bacteriocins\search\alignment_GL002679_vs_hiracinJM79.txt',
    ],
    'TableS3_public_panel': [
        W + r'\06_public_panel\metadata\panel_metadata.tsv',
        W + r'\06_public_panel\metadata\panel_metadata_all.tsv',
        W + r'\06_public_panel\metadata\download_log.tsv',
        W + r'\06_public_panel\metadata\download_record.tsv',
    ],
    'TableS4_homology_thresholds': [
        W + r'\06_public_panel\analysis\pl1_distribution.tsv',
        W + r'\06_public_panel\analysis\gene_presence_matrix.tsv',
        W + r'\06_public_panel\analysis\nisin_panel_scan.tsv',
        W + r'\06_public_panel\analysis\nisin_panel_status.tsv',
        W + r'\08_statistics\gene_frequency.tsv',
        W + r'\08_statistics\cooccurrence_fisher.tsv',
        W + r'\08_statistics\linkage_641_679.tsv',
        W + r'\08_statistics\plasmid_gene_association.tsv',
    ],
    'TableS5_MGE_and_plasmidome': [
        W + r'\04_plasmidome\annotation\mge_inventory.tsv',
        W + r'\04_plasmidome\annotation\prophage_genes.tsv',
        W + r'\04_plasmidome\annotation\prophage_annotation_full.tsv',
        W + r'\04_plasmidome\annotation\pl2_genes.tsv',
        W + r'\04_plasmidome\annotation\pl2_annotation_full.tsv',
        W + r'\04_plasmidome\annotation\pl1_key_genes.tsv',
        W + r'\05_bacteriocins\search\mge_distance_table.tsv',
    ],
    'TableS6_claim_strength': [
        W + r'\03_replicon_compare\evidence_strength.tsv',
        W + r'\02_taxonomy\taxonomy_summary.tsv',
        W + r'\03_replicon_compare\chr_difference_summary.tsv',
        W + r'\03_replicon_compare\pl2_structure_summary.tsv',
    ],
    'TableS7_functional_annotation': [
        W + r'\12_functional_annotation\summary\kegg_level1.tsv',
        W + r'\12_functional_annotation\summary\kegg_level2.tsv',
        W + r'\12_functional_annotation\summary\cog.tsv',
        W + r'\12_functional_annotation\summary\go_level2.tsv',
        W + r'\12_functional_annotation\summary\cazy_5class.tsv',
        W + r'\12_functional_annotation\summary\cazy_families.tsv',
        W + r'\12_functional_annotation\summary\card_hits.tsv',
        W + r'\12_functional_annotation\summary\ardb_hits.tsv',
        W + r'\12_functional_annotation\summary\vfdb_category_counts.tsv',
        W + r'\12_functional_annotation\summary\ani_matrix.tsv',
        W + r'\12_functional_annotation\summary\pangenome_summary.tsv',
    ],
}


def read_tsv(path):
    with open(path, encoding='utf-8', errors='replace') as f:
        return [line.rstrip('\n').split('\t') for line in f]


wb = Workbook()
wb.remove(wb.active)
for sheet, files in TABLES.items():
    ws = wb.create_sheet(sheet[:31])
    row = 1
    for path in files:
        name = os.path.basename(path)
        ws.cell(row=row, column=1, value='## ' + name).font = None
        row += 1
        if not os.path.exists(path):
            ws.cell(row=row, column=1, value='(missing: %s)' % path)
            row += 2
            continue
        if path.endswith('.txt'):
            for line in open(path, encoding='utf-8', errors='replace'):
                ws.cell(row=row, column=1, value=line.rstrip('\n'))
                row += 1
        else:
            for r_ in read_tsv(path):
                for c, v in enumerate(r_):
                    ws.cell(row=row, column=c + 1, value=v)
                row += 1
        row += 1
    # 列宽
    for col in range(1, min(ws.max_column, 12) + 1):
        ws.column_dimensions[get_column_letter(col)].width = 20

xlsx = os.path.join(OUT, 'Supplementary_Tables.xlsx')
wb.save(xlsx)
print('saved', xlsx)

# tsv 副本
for sheet, files in TABLES.items():
    dst = os.path.join(OUT, sheet)
    os.makedirs(dst, exist_ok=True)
    for path in files:
        if os.path.exists(path):
            shutil.copy2(path, os.path.join(dst, os.path.basename(path)))
print('tsv copies copied')
