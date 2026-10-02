# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
build_refs.py — 模块 05：构建 nisin 与细菌素参考集（全部记录 accession 与获取日期）
1) 11 个 nisin 系统蛋白（沿用 smbu06 分析工作区已验证 accession 集；nisR 使用 Swiss-Prot Q07597.1，
   弃用不可靠的 Q0GU37——见前期勘误记录）
2) nisin 阳性对照簇区段：L. lactis F44 CP024954.1:591,800-606,200（14,401 bp）
3) 细菌素/类IIa 参考集：UniProt 按名称检索 22 组关键词，去重
4) 阴性对照基因组由 datasets 另行下载（脚本内命令行给出）
输出：05_bacteriocins/refs/*.fasta + refs_manifest.tsv
"""
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
from p3lib import WORK, open_log

OUT = os.path.join(WORK, '05_bacteriocins', 'refs')
os.makedirs(OUT, exist_ok=True)
LOG = open_log(os.path.join(WORK, '05_bacteriocins', 'logs', 'build_refs.log'))

NISIN_ACCESSIONS = {
    'nisA': 'V5NV19', 'nisB': 'Q48673', 'nisC': 'Q48670', 'nisI': 'Q48671',
    'nisP': 'Q48674', 'nisR': 'Q07597', 'nisT': 'Q48669', 'nisF': 'Q48597',
    'nisE': 'Q48598', 'nisG': 'Q48599', 'nisK': 'Q48675',
}

BACTERIOCIN_QUERIES = [
    'enterocin', 'hiracin', 'pediocin', 'leucocin', 'sakacin', 'bavaricin',
    'piscicolin', 'divercin', 'carnobacteriocin', 'lacticin', 'lactococcin',
    'nisin', 'plantaricin', 'sublancin', 'enterolysin', 'cytolysin', 'durancin',
    'mundticin', 'garvicin', 'ubericin', 'avicin', 'bacteriocin',
]


def curl(url, out, check_fasta=True):
    r = subprocess.run(['curl', '-sL', '--max-time', '120', url], capture_output=True)
    if r.returncode != 0 or not r.stdout:
        LOG.write('FETCH FAILED: %s\n' % url)
        return False
    data = r.stdout.decode('utf-8', errors='replace')
    if check_fasta and not data.startswith('>'):
        LOG.write('NOT FASTA from %s: %s\n' % (url, data[:200]))
        return False
    with open(out, 'w') as f:
        f.write(data)
    return True


# 1) nisin 11 蛋白（逐条获取，批量 URL 会因格式报错）
p = os.path.join(OUT, 'nisin_reference_proteins.faa')
n = 0
out_lines = []
for name, acc in NISIN_ACCESSIONS.items():
    r = subprocess.run(['curl', '-sL', '--max-time', '60',
                        'https://rest.uniprot.org/uniprotkb/%s.fasta' % acc], capture_output=True)
    txt = r.stdout.decode('utf-8', errors='replace')
    if not txt.startswith('>'):
        LOG.write('nisin fetch failed: %s (%s): %s\n' % (name, acc, txt[:120]))
        continue
    head, *seq = txt.split('\n', 1)
    s = ''.join(l.strip() for l in (seq[0] if seq else '').split('\n') if l.strip() and not l.startswith('>'))
    # 取第一条记录
    s = ''.join(l.strip() for l in txt.split('\n')[1:] if l.strip() and not l.startswith('>'))
    out_lines.append('>%s|%s %s\n%s\n' % (name, acc, head.split(' ', 1)[1][:70] if ' ' in head else '', s))
    n += 1
with open(p, 'w') as f:
    f.write(''.join(out_lines))
LOG.write('nisin proteins fetched: %d (expect 11)\n' % n)

# 2) 阳性对照簇区段
p2 = os.path.join(OUT, 'nisin_positive_cluster_F44.fasta')
ok2 = curl('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP024954.1'
           '&seq_start=591800&seq_stop=606200&rettype=fasta&retmode=text', p2)
LOG.write('F44 cluster region fetched: %s\n' % ok2)

# 3) 细菌素参考集（search 端点 + 长度 ≤400 aa 过滤，控噪）
p3 = os.path.join(OUT, 'bacteriocin_refset.faa')
seen = set()
seqs = []
for q in BACTERIOCIN_QUERIES:
    url = ('https://rest.uniprot.org/uniprotkb/search?query=%s+AND+length:%%5B1+TO+400%%5D'
           '&format=fasta&size=500' % q)
    r = subprocess.run(['curl', '-sL', '--max-time', '60', url], capture_output=True)
    if r.returncode != 0 or not r.stdout:
        LOG.write('query failed: %s\n' % q)
        continue
    txt = r.stdout.decode('utf-8', errors='replace')
    cnt = 0
    for block in txt.split('>')[1:]:
        head, *seq = block.split('\n')
        s = ''.join(l.strip() for l in seq if l.strip())
        if s and s not in seen:
            seen.add(s)
            acc = head.split('|')[1] if '|' in head else head.split()[0]
            seqs.append('>%s %s\n%s\n' % (acc, head.split(' ', 1)[1][:80] if ' ' in head else '', s))
            cnt += 1
    LOG.write('query %s: +%d\n' % (q, cnt))
with open(p3, 'w') as f:
    f.write(''.join(seqs))
LOG.write('bacteriocin refset total: %d unique sequences\n' % len(seqs))

# refs manifest
import datetime
with open(os.path.join(WORK, '05_bacteriocins', 'refs_manifest.tsv'), 'w', encoding='utf-8') as f:
    f.write('file\tsource\taccessions_or_query\tfetch_date\tn_sequences\n')
    f.write('nisin_reference_proteins.faa\tUniProt REST\t%s\t%s\t%d\n' % (
        ';'.join('%s=%s' % (k, v) for k, v in NISIN_ACCESSIONS.items()),
        datetime.date.today().isoformat(), n))
    f.write('nisin_positive_cluster_F44.fasta\tNCBI nuccore CP024954.1:591800-606200\tCP024954.1\t%s\t1\n' % datetime.date.today().isoformat())
    f.write('bacteriocin_refset.faa\tUniProt REST stream\t%s\t%s\t%d\n' % (
        ';'.join(BACTERIOCIN_QUERIES), datetime.date.today().isoformat(), len(seqs)))
print('build_refs done: nisin=%d, refset=%d' % (n, len(seqs)))
