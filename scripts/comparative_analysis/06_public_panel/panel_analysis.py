# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
panel_analysis.py — 模块 06：公共面板元数据、去重、共享质粒分布与候选基因存在/缺失
A) 元数据：从各下载包 assembly_data_report.jsonl 重建完整 panel_metadata_all.tsv
B) 去重：GCA/GCF 同 assembly 重复 + (organism, strain) 重复保留最优 → panel_metadata.tsv
C) 分布：pl1 全长（blastn）分层统计 + 模块区分布
D) 基因存在/缺失：tblastn 候选 3 蛋白 + nisin 6 特异蛋白；主阈值 id≥70/qcov≥70（另 60/80 敏感性）
输出：06_public_panel/analysis/*.tsv
"""
import glob
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
import numpy as np
from p3lib import (WORK, IN06, IN08, BLASTN, TBLASTN, MAKEDB, DATASETS, read_fasta, write_fasta,
                   revcomp, run, open_log)

OUT = os.path.join(WORK, '06_public_panel')
AN = os.path.join(OUT, 'analysis')
os.makedirs(AN, exist_ok=True)
LOG = open_log(os.path.join(OUT, 'logs', 'panel_analysis.log'))

# ---------- A) 元数据（datasets summary 批量；jsonl 作后备） ----------
acc_list = [l.strip() for l in open(os.path.join(OUT, 'metadata', 'panel_accessions.txt'), encoding='utf-8') if l.strip()]
meta = []
CH = 40
for i in range(0, len(acc_list), CH):
    chunk = acc_list[i:i + CH]
    r = subprocess.run([DATASETS, 'summary', 'genome', 'accession'] + chunk + ['--as-json-lines'],
                       capture_output=True, text=True, encoding='utf-8', errors='replace', timeout=600)
    for line in (r.stdout or '').splitlines():
        line = line.strip()
        if not line.startswith('{'):
            continue
        d = json.loads(line)
        ai = d.get('assembly_info', {})
        bs = ai.get('biosample') or {}
        o = d.get('organism', {}) or {}
        st = d.get('assembly_stats', {}) or {}
        meta.append(dict(
            accession=d.get('accession', ''),
            organism=o.get('organism_name', ''),
            strain=(o.get('infraspecific_names') or {}).get('strain', '') or bs.get('strain', ''),
            assembly_name=ai.get('assembly_name', ''),
            assembly_level=ai.get('assembly_level', ''),
            total_length_bp=st.get('total_sequence_length', ''),
            gc_percent=st.get('gc_percent', ''),
            contig_n50=st.get('contig_n50', ''),
            number_of_contigs=st.get('number_of_contigs', ''),
            isolation_source=bs.get('isolation_source', '') or '',
            host=(o.get('infraspecific_names') or {}).get('host', '') or '',
            geo_loc_name=bs.get('geo_loc_name', '') or '',
            collection_date=bs.get('collection_date', '') or '',
            biosample=bs.get('accession', ''),
            bioproject=ai.get('bioproject_accession', ''),
            assembly_status=ai.get('assembly_status', ''),
            annotated_protein_coding=(d.get('annotation_info', {}) or {}).get('stats', {}).get('gene_counts', {}).get('protein_coding', ''),
            submitter=ai.get('submitter', ''),
        ))
    LOG.write('summary chunk %d: %d records\n' % (i // CH, len(meta)))
LOG.write('metadata entries: %d (from %d accessions)\n' % (len(meta), len(acc_list)))

# 去重：same assembly_name (GCA/GCF 重复) + same (organism, strain) 保留最优
by_key = {}
rank = {'Complete Genome': 0, 'Chromosome': 1, 'Scaffold': 2, 'Contig': 3}
for m in meta:
    key = (m['organism'].lower(), m['strain'].lower() or m['assembly_name'])
    m['_rank'] = rank.get(m['assembly_level'], 9)
    cur = by_key.get(key)
    if cur is None or (m['_rank'], -int(m['total_length_bp'] or 0)) < (cur['_rank'], -int(cur['total_length_bp'] or 0)):
        by_key[key] = m
dedup = list(by_key.values())
LOG.write('dedup: %d -> %d\n' % (len(meta), len(dedup)))

cols = ['accession', 'organism', 'strain', 'assembly_name', 'assembly_level', 'total_length_bp',
        'gc_percent', 'contig_n50', 'number_of_contigs', 'isolation_source', 'host',
        'geo_loc_name', 'collection_date', 'biosample', 'bioproject', 'assembly_status',
        'annotated_protein_coding', 'submitter']
for path, rows_ in [(os.path.join(OUT, 'metadata', 'panel_metadata_all.tsv'), meta),
                    (os.path.join(OUT, 'metadata', 'panel_metadata.tsv'), dedup)]:
    with open(path, 'w', encoding='utf-8') as f:
        f.write('\t'.join(cols) + '\n')
        for m in sorted(rows_, key=lambda x: x['accession']):
            f.write('\t'.join(str(m.get(c, '')) for c in cols) + '\n')

# 下载记录（date, URL）
import datetime
with open(os.path.join(OUT, 'metadata', 'download_record.tsv'), 'w', encoding='utf-8') as f:
    f.write('accession\tdownload_date\tsource\n')
    today = datetime.date.today().isoformat()
    for m in sorted(meta, key=lambda x: x['accession']):
        f.write('%s\t%s\thttps://www.ncbi.nlm.nih.gov/datasets/genome/%s/\n' % (m['accession'], today, m['accession']))

# ---------- C/D) 组合 DB ----------
combined = os.path.join(AN, 'panel_all.fna')
n = 0
CONTIG_LEN = {}
with open(combined, 'w') as fo:
    for m in sorted(dedup, key=lambda x: x['accession']):
        d = os.path.join(OUT, 'downloads', m['accession'])
        fna = None
        for root, _, files in os.walk(d):
            for fn in files:
                if fn.endswith('.fna'):
                    fna = os.path.join(root, fn)
                    break
            if fna:
                break
        if not fna:
            LOG.write('missing fna: %s\n' % m['accession'])
            continue
        for name, s in read_fasta(fna):
            key = '%s|%s' % (m['accession'], name.split()[0])
            fo.write('>%s\n%s\n' % (key, s))
            CONTIG_LEN[key] = len(s)
            n += 1
LOG.write('combined: %d contigs from %d genomes\n' % (n, len(dedup)))
db = os.path.join(AN, 'panel_db')
if not os.path.exists(db + '.nsq'):
    run([MAKEDB, '-in', combined, '-dbtype', 'nucl', '-out', db], log=LOG, check=True)
else:
    LOG.write('reuse existing panel_db\n')

# pl1 全长分布
pl1 = read_fasta(WORK + r'\03_replicon_compare\alignments\smbu06_pl1.fa')[0][1]
pl1fa = os.path.join(AN, 'pl1_query.fa')
write_fasta(pl1fa, [('smbu06_Plasmid1', pl1)])
r = run([BLASTN, '-task', 'megablast', '-query', pl1fa, '-db', db,
         '-outfmt', '6 sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen',
         '-evalue', '1e-5', '-max_target_seqs', '5000', '-max_hsps', '100', '-num_threads', '8'],
        log=LOG, check=True)
from collections import defaultdict
# 13 列: sseqid(0) pident(1) length(2) mismatch(3) gapopen(4) qstart(5) qend(6) sstart(7) send(8) ...
raw = defaultdict(list)
for line in r.stdout.strip().split('\n'):
    if not line:
        continue
    p = line.split('\t')
    acc = p[0].split('|')[0]
    raw[acc].append(dict(subject=p[0], pid=float(p[1]), alen=int(p[2]),
                         qs=int(p[5]), qe=int(p[6])))
per_g = {}
for acc, hits in raw.items():
    # 全 contig 汇总，query 坐标共线链（跨 contig 去重；被链覆盖的内部重复命中剔除）
    hs = sorted(hits, key=lambda h: (min(h['qs'], h['qe']), -h['alen']))
    d = dict(cov=np.zeros(len(pl1), dtype=bool), idents=[], alen=0, contigs=set(), main_contig='', main_len=0)
    end = 0
    for h in hs:
        a, b = min(h['qs'], h['qe']) - 1, max(h['qs'], h['qe'])
        if b <= end:
            continue
        d['cov'][a:b] = True
        d['idents'].append((h['pid'], h['alen']))
        d['alen'] += h['alen']
        d['contigs'].add(h['subject'])
        if h['alen'] > d['main_len']:
            d['main_len'], d['main_contig'] = h['alen'], h['subject']
        end = max(end, b)
    per_g[acc] = d
rows = []
for acc, d in sorted(per_g.items()):
    covp = d['cov'].sum() / len(pl1) * 100
    wt = sum(i * l for i, l in d['idents']) / max(1, sum(l for _, l in d['idents']))
    if covp >= 95 and wt >= 99:
        cls = 'full_length_near_identical'
    elif covp >= 95:
        cls = 'full_length_divergent'
    elif covp >= 30:
        cls = 'partial_backbone'
    elif covp >= 5:
        cls = 'fragment'
    else:
        cls = 'trace'
    rows.append((acc, round(covp, 1), round(wt, 2), len(d['contigs']), cls,
                 '%s(%dbp)' % (d['main_contig'], CONTIG_LEN.get(d['main_contig'], 0))))

with open(os.path.join(AN, 'pl1_distribution.tsv'), 'w', encoding='utf-8') as f:
    f.write('accession\tquery_cov_pct\tidentity_wt\tcontigs\tclass\tmain_contig\n')
    for row in rows:
        f.write('\t'.join(str(x) for x in row) + '\n')
LOG.write('pl1 distribution: %d genomes with hits; classes: %s\n' % (
    len(rows), {c: sum(1 for r_ in rows if r_[4] == c) for c in set(r_[4] for r_ in rows)}))

# 模块区（641+642+679 及邻域 31–75 kb）分布
module_seq = pl1[31073 - 1:74688]
write_fasta(os.path.join(AN, 'module_region.fa'), [('module_31-75kb', module_seq)])

# 基因存在/缺失（tblastn）
GENES = {
    'enterocinP_like_GL002641': os.path.join(WORK, '05_bacteriocins', 'search', 'smbu06_proteins.faa'),
    'hiracinJM79_like_GL002679': None,
    'sakacinA_immunity_GL002642': None,
}
prot = dict((x.split()[0], s) for x, s in read_fasta(os.path.join(WORK, '05_bacteriocins', 'search', 'smbu06_proteins.faa')))
qfa = os.path.join(AN, 'candidates.faa')
write_fasta(qfa, [('enterocinP_like', prot['smbu06GL002641']),
                  ('hiracinJM79_like', prot['smbu06GL002679']),
                  ('sakacinA_immunity', prot['smbu06GL002642'])])
r = run([TBLASTN, '-query', qfa, '-db', db, '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen',
         '-evalue', '1e-5', '-seg', 'no', '-max_target_seqs', '5000', '-max_hsps', '3', '-num_threads', '8'],
        log=LOG, check=True)
hits = defaultdict(list)
with open(os.path.join(AN, 'gene_hits_all.tsv'), 'w', encoding='utf-8') as fh:
    fh.write('gene\taccession\tcontig\tpident\talen\tqcov\tevalue\tsstart\tsend\n')
    for line in r.stdout.strip().split('\n'):
        if not line:
            continue
        p = line.split('\t')
        acc = p[1].split('|')[0]
        qcov = (int(p[7]) - int(p[6]) + 1) / float(p[12]) * 100
        hits[(p[0], acc)].append(dict(pid=float(p[2]), qcov=qcov, evalue=float(p[10])))
        fh.write('%s\t%s\t%s\t%s\t%s\t%.1f\t%s\t%s\t%s\n' % (
            p[0].split('|')[0], acc, p[1], p[2], p[3], qcov, p[10], p[8], p[9]))

with open(os.path.join(AN, 'gene_presence_matrix.tsv'), 'w', encoding='utf-8') as f:
    f.write('accession\tgene\tbest_pident\tbest_qcov\tpresence_70_70\tpresence_60_70\tpresence_80_80\n')
    for m in sorted(dedup, key=lambda x: x['accession']):
        for gene in ['enterocinP_like', 'hiracinJM79_like', 'sakacinA_immunity']:
            hs = hits.get((gene, m['accession']), [])
            best = max(hs, key=lambda h: (h['pid'] * h['qcov'])) if hs else None
            f.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
                m['accession'], gene,
                round(best['pid'], 1) if best else '',
                round(best['qcov'], 1) if best else '',
                'yes' if best and best['pid'] >= 70 and best['qcov'] >= 70 else 'no',
                'yes' if best and best['pid'] >= 60 and best['qcov'] >= 70 else 'no',
                'yes' if best and best['pid'] >= 80 and best['qcov'] >= 80 else 'no'))
LOG.write('gene presence matrix written\n')

# nisin 面板扫描（6 特异蛋白）
nis_faa = os.path.join(WORK, '05_bacteriocins', 'refs', 'nisin_reference_proteins.faa')
with open(os.path.join(AN, 'nisin_panel_scan.tsv'), 'w', encoding='utf-8') as f:
    f.write('accession\tgene\tpident\tqcov\tevalue\tstatus\n')
    r = run([TBLASTN, '-query', nis_faa, '-db', db,
             '-outfmt', '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen',
             '-evalue', '1e-5', '-seg', 'no', '-max_target_seqs', '600', '-max_hsps', '1', '-num_threads', '8'],
            log=LOG, check=True)
    best = {}
    for line in r.stdout.strip().split('\n'):
        if not line:
            continue
        p = line.split('\t')
        g = p[0].split('|')[0]
        acc = p[1].split('|')[0]
        qcov = (int(p[5]) - int(p[4]) + 1) / float(p[10]) * 100
        k = (g, acc)
        if k not in best or float(p[2]) > best[k][0]:
            best[k] = (float(p[2]), qcov, float(p[10]))
    hit_genes = {}
    for (g, acc), (pid, qcov, ev) in sorted(best.items()):
        st = 'hit' if pid >= 80 and qcov >= 80 else ('weak' if pid >= 40 else 'noise')
        f.write('%s\t%s\t%.1f\t%.1f\t%.1e\t%s\n' % (acc, g, pid, qcov, ev, st))
        if st == 'hit':
            hit_genes.setdefault(acc, set()).add(g)
    # 汇总
    with open(os.path.join(AN, 'nisin_panel_status.tsv'), 'w', encoding='utf-8') as f2:
        f2.write('accession\tnisin_specific_hits\tcluster_call\n')
        for m in sorted(dedup, key=lambda x: x['accession']):
            spec = sorted(hit_genes.get(m['accession'], set()) & {'nisA', 'nisB', 'nisC', 'nisI', 'nisP', 'nisR'})
            call = 'complete_cluster' if all(g in spec for g in ('nisA', 'nisB', 'nisC')) else ('partial' if spec else 'not_detected')
            f2.write('%s\t%s\t%s\n' % (m['accession'], ';'.join(spec), call))
LOG.write('nisin panel scan written\n')
print('panel_analysis done')
