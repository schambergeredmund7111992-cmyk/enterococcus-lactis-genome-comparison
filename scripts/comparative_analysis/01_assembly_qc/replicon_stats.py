# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
replicon_stats.py — 模块 01：复制子级组装统计（Table S1 基础）
来源:交付 FASTA(长度/GC/N) + 公司表(topology/覆盖/深度/基因/ncRNA/校正) + GFF(CDS 计数)
输出:01_assembly_qc/replicon_stats.tsv
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
from p3lib import WORK, IN06, IN08, read_fasta, gc_content, open_log

OUT = os.path.join(WORK, '01_assembly_qc', 'replicon_stats.tsv')
LOG = open_log(os.path.join(WORK, '01_assembly_qc', 'logs', 'replicon_stats.log'))

STRAINS = {
    'smbu06': {
        'root': IN06 + r'\smbu06',
        'fastas': {'Chromosome1': IN06 + r'\smbu06\2.Assembly\smbu06.Chromosome.fasta',
                   'Plasmid1': IN06 + r'\smbu06\2.Assembly\smbu06.Plasmid.fasta'},
    },
    'smbu08': {
        'root': IN08 + r'\smbu08',
        'fastas': {'Chromosome1': IN08 + r'\smbu08\2.Assembly\smbu08.Chromosome.fasta',
                   'Plasmid1': None,  # 从合并文件拆
                   'Plasmid2': None},
    },
}


def tsv_rows(path):
    rows = []
    with open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.rstrip('\n').rstrip('\r')
            if line.strip() == '':
                continue
            rows.append(line.split('\t'))
    return rows


def num(s):
    return s.replace(',', '').strip()


def parse_topology(root, prefix):
    out = {}
    for r in tsv_rows(os.path.join(root, '2.Assembly', prefix + '.Cir.Stat.xls'))[1:]:
        if len(r) >= 2:
            out[r[0]] = 'circular' + (' (speculative)' if len(r) > 2 and r[2].strip() else '')
    return out


def parse_coverage(root, prefix):
    out = {}
    rows = tsv_rows(os.path.join(root, '2.Assembly', prefix + '.CoverageStatCom.xls'))
    for r in rows[1:]:
        if len(r) >= 7 and r[0].strip() or (len(r) >= 7 and r[1].strip()):
            name = r[1].strip()
            if r[0].strip():
                out['__sample_depth__'] = r[5].strip()
                out['__usage__'] = r[6].strip()
            if name and name != 'Total':
                out[name] = {'coverage_pct': num(r[4]), 'depth': num(r[5])}
    return out


def parse_gene_stat(root, prefix):
    rows = tsv_rows(os.path.join(root, '3.Genome_Component', 'Gene_Predict', prefix + '.Gene.stat.xls'))
    r = rows[1]
    return {'genes_total': num(r[2]), 'gene_len': num(r[3]), 'gene_gc': r[6].strip()}


def parse_ncrna(root, prefix):
    out = {}
    rows = tsv_rows(os.path.join(root, '3.Genome_Component', 'ncRNA_Finding', prefix + '.ncRNA.stat.xls'))
    for r in rows[1:]:
        if len(r) >= 3:
            out[r[1].strip()] = num(r[2])
    return out


def parse_correct(root, prefix):
    rows = tsv_rows(os.path.join(root, '2.Assembly', prefix + '.CorrectRate.stat.xls'))
    r = rows[1]
    return {'indel': num(r[2]), 'snp': num(r[3]), 'rate_pct': r[5].strip()}


def parse_kmer(root, prefix):
    rows = tsv_rows(os.path.join(root, '2.Assembly', prefix + '.kmer.stat.xls'))
    r = rows[1]
    return {'k': r[1], 'pk_depth': num(r[3]), 'kmer_genome_mb': num(r[4])}


def count_cds_per_replicon(root, prefix):
    counts = {}
    path = os.path.join(root, '3.Genome_Component', 'Gene_Predict', prefix + '.Gene.gff')
    with open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            p = line.split('\t')
            if len(p) > 2 and p[2] == 'CDS':
                counts[p[0]] = counts.get(p[0], 0) + 1
    return counts


def fstats(seq):
    n = len(seq)
    n_count = seq.count('N') + seq.count('n')
    return n, gc_content(seq), n_count


def main():
    out_rows = []
    for strain, cfg in STRAINS.items():
        root, prefix = cfg['root'], strain
        topo = parse_topology(root, prefix)
        cov = parse_coverage(root, prefix)
        gs = parse_gene_stat(root, prefix)
        ncrna = parse_ncrna(root, prefix)
        corr = parse_correct(root, prefix)
        km = parse_kmer(root, prefix)
        cds_counts = count_cds_per_replicon(root, prefix)
        LOG.write('%s: topo=%s cov=%s\n' % (strain, topo, cov))

        # 复制子序列
        seqs = {}
        if strain == 'smbu06':
            for name, fa in cfg['fastas'].items():
                recs = read_fasta(fa)
                seqs[name] = recs[0][1]
        else:
            recs = read_fasta(cfg['root'] + r'\2.Assembly\smbu08.Chromosome.fasta')
            seqs['Chromosome1'] = recs[0][1]
            for r in read_fasta(IN08 + r'\smbu08\2.Assembly\smbu08.Plasmid.fasta'):
                if 'Plasmid1' in r[0]:
                    seqs['Plasmid1'] = r[1]
                elif 'Plasmid2' in r[0]:
                    seqs['Plasmid2'] = r[1]

        for name in ['Chromosome1', 'Plasmid1', 'Plasmid2']:
            if name not in seqs:
                continue
            n, gc, ncnt = fstats(seqs[name])
            c = cov.get(name, {})
            out_rows.append([
                strain, name, n, '%.2f' % gc, ncnt,
                topo.get(name, ''),
                c.get('coverage_pct', ''), c.get('depth', ''),
                cds_counts.get(name, ''),
                ncrna.get('tRNA', '') if name == 'Chromosome1' else '',
                (ncrna.get('16s_rRNA (Denovo)', ''), ncrna.get('23s_rRNA (Denovo)', ''), ncrna.get('5s_rRNA (Denovo)', ''))
                if name == 'Chromosome1' else '',
                ncrna.get('sRNA', '') if name == 'Chromosome1' else '',
                corr['indel'] if name == 'Chromosome1' else '',
                corr['snp'] if name == 'Chromosome1' else '',
                corr['rate_pct'] if name == 'Chromosome1' else '',
                km['pk_depth'] if name == 'Chromosome1' else '',
            ])
        # reproris overall
        LOG.write('%s genes_total=%s len=%s ncrna=%s\n' % (strain, gs, gs.get('gene_len'), ncrna))

    with open(OUT, 'w', encoding='utf-8') as f:
        f.write('strain\treplicon\tlength_bp\tgc_pct\tn_bases\ttopology\tcoverage_pct_reads(delivered)\tdepth(delivered)\tcds_count(gff)\ttrna_genes\trrna_16s_23s_5s\tsrna\tindel_corrected\tsnp_corrected\tcorrect_rate_pct\tkmer_peak_depth\n')
        for r in out_rows:
            f.write('\t'.join(str(x) for x in r) + '\n')
    print('wrote', OUT)
    LOG.write('DONE\n')


if __name__ == '__main__':
    main()
