# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
taxonomy_check.py — 模块 02：计算分类核验（只作核验，不修改项目命名）
1) 参考基因组消息与元数据记录
2) 字节级环状比较：smbu06/smbu08 各复制子 vs E. lactis 194 (GCF_056582645.1)
3) ANIb（JSpeciesWS 口径 1020bp 片段，双向）vs 194 / L. lactis 14B4 / L. lactis MG1363
4) 16S：本组装 16S 拷贝 vs 194 及外部 EF102815(EF102815 声称来源) / EF100778(F-119 真实 16S)
输出：02_taxonomy/taxonomy_summary.tsv、ani_results.tsv、16S_comparison.tsv、ref_metadata.json
"""
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '00_manifest', 'scripts'))
import numpy as np
from p3lib import (WORK, IN06, IN08, BLASTN, MAKEDB, read_fasta, write_fasta,
                   revcomp, sha256_text, canonical_circular, run, open_log)

OUT = os.path.join(WORK, '02_taxonomy')
REF = os.path.join(OUT, 'refs')
TMP = os.path.join(OUT, 'tmp')
os.makedirs(TMP, exist_ok=True)
LOG = open_log(os.path.join(OUT, 'logs', 'taxonomy_check.log'))

REF_FILES = {
    'E_lactis_194': {
        'chr': REF + r'\GCF_056582645.1\ncbi_dataset\data\GCF_056582645.1\GCF_056582645.1_ASM5658264v1_genomic.fna',
        'json': REF + r'\GCF_056582645.1\ncbi_dataset\data\assembly_data_report.jsonl',
    },
    'L_lactis_14B4': {
        'chr': REF + r'\GCF_003176835.1\ncbi_dataset\data\GCF_003176835.1\GCF_003176835.1_ASM317683v1_genomic.fna',
    },
    'L_lactis_MG1363': {
        'chr': REF + r'\GCF_000009425.1\ncbi_dataset\data\GCF_000009425.1\GCF_000009425.1_ASM942v1_genomic.fna',
    },
}

# ---------- 载入 ----------
refs = {}   # name -> {'chromosome': seq, 'plasmid': seq}
for name, cfg in REF_FILES.items():
    recs = read_fasta(cfg['chr'])
    refs[name] = {}
    for n, s in recs:
        key = 'chromosome' if 'chromosome' in n.lower() else ('plasmid' if 'plasmid' in n.lower() else n)
        refs[name][key] = s
    LOG.write('%s: %s\n' % (name, {k: len(v) for k, v in refs[name].items()}))

asm = {
    'smbu06_chr': read_fasta(IN06 + r'\smbu06\2.Assembly\smbu06.Chromosome.fasta')[0][1],
    'smbu06_pl1': read_fasta(IN06 + r'\smbu06\2.Assembly\smbu06.Plasmid.fasta')[0][1],
    'smbu08_chr': read_fasta(IN08 + r'\smbu08\2.Assembly\smbu08.Chromosome.fasta')[0][1],
}
pl08 = read_fasta(IN08 + r'\smbu08\2.Assembly\smbu08.Plasmid.fasta')
asm['smbu08_pl1'] = [r for r in pl08 if 'Plasmid1' in r[0]][0][1]

# 元数据
meta = {}
with open(REF_FILES['E_lactis_194']['json'], encoding='utf-8') as f:
    for line in f:
        d = json.loads(line)
        a = d.get('assemblyInfo', {})
        o = d.get('organism', {})
        meta = {
            'organism': o.get('organismName'),
            'infraspecific': o.get('infraspecificNames', {}),
            'assembly_name': a.get('assemblyName'),
            'assembly_level': a.get('assemblyLevel'),
            'bioproject': a.get('bioprojectAccession'),
            'biosample': a.get('biosampleAccession'),
            'submitter': a.get('submitter'),
            'submission_date': a.get('submissionDate'),
            'refseq_category': a.get('refseqCategory'),
            'contig_n50': a.get('contigN50'),
        }
        break
with open(os.path.join(OUT, 'ref_metadata.json'), 'w', encoding='utf-8') as f:
    json.dump(meta, f, ensure_ascii=False, indent=1)
LOG.write('194 metadata: %s\n' % meta)

# ---------- 字节级环状比较 ----------
def circular_relation(A, B):
    """返回 ('identical'|'identical_rc'|'other', offset, note)"""
    if A == B:
        return 'identical', 0, 'byte-identical linear representation'
    t = A + A
    r = t.find(B)
    if r >= 0:
        return 'identical', r, 'forward rotation by %d bp' % r
    A2 = A + A
    Brc = revcomp(B)
    r2 = A2.find(Brc)
    if r2 >= 0:
        return 'identical_rc', r2, 'reverse-complement rotation by %d bp' % r2
    return 'other', None, 'not identical as circular molecules'


# 用 blastn 定位并统计差异（near-identical 情形）
def block_diff(qname, Q, sname, S, tag):
    qfa = os.path.join(TMP, '%s.fa' % qname)
    write_fasta(qfa, [(qname, Q)])
    sfa = os.path.join(TMP, '%s.fa' % sname)
    write_fasta(sfa, [(sname, S)])
    db = os.path.join(TMP, '%s_db' % sname)
    run([MAKEDB, '-in', sfa, '-dbtype', 'nucl', '-out', db], log=LOG, check=True)
    r = run([BLASTN, '-task', 'megablast', '-query', qfa, '-db', db,
             '-outfmt', '6 pident length mismatch gapopen qstart qend sstart send bitscore',
             '-evalue', '1e-5', '-max_hsps', '2000', '-num_threads', '8'], log=LOG, check=True)
    hsps = []
    for line in r.stdout.strip().split('\n'):
        if not line.strip():
            continue
        p = line.split('\t')
        hsps.append(dict(pident=float(p[0]), alen=int(p[1]), mm=int(p[2]), gapopen=int(p[3]),
                         qs=int(p[4]), qe=int(p[5]), ss=int(p[6]), se=int(p[7]), bits=float(p[8])))
    big = sorted([h for h in hsps if h['alen'] >= 50000 and h['ss'] < h['se']], key=lambda h: h['qs'])
    # 差异统计
    events = []
    cov = np.zeros(len(Q), dtype=bool)
    for h in big:
        cov[h['qs'] - 1:h['qe']] = True
        a = Q[h['qs'] - 1:h['qe']]
        b = S[h['ss'] - 1:h['se']]
        if len(a) == len(b):
            aa = np.frombuffer(a.encode(), dtype=np.uint8)
            bb = np.frombuffer(b.encode(), dtype=np.uint8)
            mm = np.nonzero(aa != bb)[0]
            for i in mm:
                events.append(('SNP', h['qs'] + int(i), h['ss'] + int(i), a[i], b[i]))
        else:
            pre = 0
            while pre < min(len(a), len(b)) and a[pre] == b[pre]:
                pre += 1
            suf = 0
            while suf < min(len(a), len(b)) - pre and a[len(a) - 1 - suf] == b[len(b) - 1 - suf]:
                suf += 1
            events.append(('INDEL_CORE', h['qs'] + pre, h['ss'] + pre,
                           'q=%dbp' % (len(a) - pre - suf), 's=%dbp' % (len(b) - pre - suf)))
    # S 中未覆盖区
    covs = np.zeros(len(S), dtype=bool)
    for h in big:
        covs[h['ss'] - 1:h['se']] = True
    def gaps(mask):
        idx = np.nonzero(~mask)[0]
        out = []
        if len(idx):
            st = prev = idx[0]
            for x in idx[1:]:
                if x != prev + 1:
                    out.append((int(st) + 1, int(prev) + 1))
                    st = x
                prev = x
            out.append((int(st) + 1, int(prev) + 1))
        return out
    snps = sum(1 for e in events if e[0] == 'SNP')
    LOG.write('%s: big_hsps=%d q_cov=%.2f%% snps=%d cores=%d uncovered_in_subject=%s\n' % (
        tag, len(big), cov.sum() / len(Q) * 100, snps,
        sum(1 for e in events if e[0] != 'SNP'), gaps(covs)))
    return dict(big=len(big), q_cov_pct=round(cov.sum() / len(Q) * 100, 4), snps=snps,
                cores=[e for e in events if e[0] != 'SNP'], uncovered_subject=gaps(covs),
                min_pid=min(h['pident'] for h in big) if big else None,
                wt_pid=(sum(h['pident'] * h['alen'] for h in big) / sum(h['alen'] for h in big)) if big else None)

summary_rows = []
for asm_name, asm_seq, ref_name, ref_key in [
    ('smbu06_chr', asm['smbu06_chr'], 'E_lactis_194', 'chromosome'),
    ('smbu06_pl1', asm['smbu06_pl1'], 'E_lactis_194', 'plasmid'),
    ('smbu08_chr', asm['smbu08_chr'], 'E_lactis_194', 'chromosome'),
    ('smbu08_pl1', asm['smbu08_pl1'], 'E_lactis_194', 'plasmid'),
]:
    R = refs[ref_name][ref_key]
    rel, off, note = circular_relation(R, asm_seq)
    LOG.write('%s vs %s/%s: %s (%s)\n' % (asm_name, ref_name, ref_key, rel, note))
    row = dict(assembly=asm_name, reference='%s/%s' % (ref_name, ref_key),
               len_asm=len(asm_seq), len_ref=len(R), relation=rel, note=note)
    if rel == 'other' and len(asm_seq) > 100000 and abs(len(asm_seq) - len(R)) < 100000:
        bd = block_diff(asm_name, asm_seq, ref_name + '_' + ref_key, R, '%s_vs_%s' % (asm_name, ref_name))
        row.update(bd)
    summary_rows.append(row)
    print(asm_name, 'vs', ref_name, ref_key, ':', rel, note)

with open(os.path.join(OUT, 'taxonomy_summary.tsv'), 'w', encoding='utf-8') as f:
    keys = ['assembly', 'reference', 'len_asm', 'len_ref', 'relation', 'note', 'big', 'q_cov_pct',
            'snps', 'min_pid', 'wt_pid', 'uncovered_subject', 'cores']
    f.write('\t'.join(keys) + '\n')
    for r in summary_rows:
        f.write('\t'.join(str(r.get(k, '')) for k in keys) + '\n')

# ---------- ANIb ----------
FRAG = 1020


def frag_write(seqs, out):
    n = 0
    with open(out, 'w') as f:
        for name, s in seqs:
            for i in range(0, len(s) - FRAG + 1, FRAG):
                n += 1
                f.write('>%s_frag%d\n%s\n' % (name, n, s[i:i + FRAG]))
    return n


def ani_one(query_frags, db):
    r = run([BLASTN, '-query', query_frags, '-db', db,
             '-outfmt', '6 qseqid pident length qlen', '-evalue', '1e-5',
             '-max_target_seqs', '5', '-num_threads', '8'], log=None)
    best = {}
    for line in r.stdout.strip().split('\n'):
        if not line:
            continue
        qid, pid, alen, qlen = line.split('\t')
        if qid not in best:
            best[qid] = (float(pid), int(alen), int(qlen))
    total = sum(1 for _ in open(query_frags)) // 2
    passed = [(p, a) for p, a, q in best.values() if p >= 70 and a >= 0.7 * q]
    if not passed:
        return None, 0, total
    tl = sum(a for _, a in passed)
    return sum(p * a for p, a in passed) / tl, len(passed), total


ani_rows = []
for asm_tag, asm_seq in [('smbu06', asm['smbu06_chr']), ('smbu08', asm['smbu08_chr'])]:
    afa = os.path.join(TMP, 'asm_%s.fa' % asm_tag)
    write_fasta(afa, [(asm_tag, asm_seq)])
    adb = os.path.join(TMP, 'db_asm_%s' % asm_tag)
    run([MAKEDB, '-in', afa, '-dbtype', 'nucl', '-out', adb], log=None, check=True)
    qf = os.path.join(TMP, 'frags_%s.fa' % asm_tag)
    nf = frag_write([(asm_tag, asm_seq)], qf)
    LOG.write('ANI frags %s: %d (%dbp each)\n' % (asm_tag, nf, FRAG))
    for ref_name, ref_d in refs.items():
        rc = ref_d.get('chromosome')
        if not rc:
            continue
        rfa = os.path.join(TMP, 'ref_%s.fa' % ref_name)
        write_fasta(rfa, [(ref_name, rc)])
        rdb = os.path.join(TMP, 'db_ref_%s' % ref_name)
        run([MAKEDB, '-in', rfa, '-dbtype', 'nucl', '-out', rdb], log=None, check=True)
        rf = os.path.join(TMP, 'frags_ref_%s.fa' % ref_name)
        frag_write([(ref_name, rc)], rf)
        a_ab, na, ta = ani_one(qf, rdb)
        a_ba, nb, tb = ani_one(rf, adb)
        vals = [v for v in (a_ab, a_ba) if v]
        ani = sum(vals) / len(vals) if vals else 0.0
        ani_rows.append((asm_tag, ref_name, round(ani, 2), '%d/%d' % (na, ta), '%d/%d' % (nb, tb)))
        print('ANI %s vs %s = %.2f%%' % (asm_tag, ref_name, ani))
with open(os.path.join(OUT, 'ani_results.tsv'), 'w', encoding='utf-8') as f:
    f.write('query\treference\tANI_pct\tquery_frags_pass\tref_frags_pass\n')
    for r in ani_rows:
        f.write('\t'.join(map(str, r)) + '\n')

# ---------- 16S ----------
# 本组装 16S 取自交付的 de novo rRNA 预测
def load_16s(path):
    recs = read_fasta(path)
    return [(n, s) for n, s in recs if '16s' in n.lower()]

s06_16s = load_16s(IN06 + r'\smbu06\3.Genome_Component\ncRNA_Finding\smbu06.denovo.rRNA.fasta')
s08_16s = load_16s(IN08 + r'\smbu08\3.Genome_Component\ncRNA_Finding\smbu08.denovo.rRNA.fasta')
LOG.write('16S copies: smbu06=%d smbu08=%d\n' % (len(s06_16s), len(s08_16s)))

# 外部 16S 序列（EF102815 声称来源；EF100778 F-119 真实 16S）
ext16 = {}
for acc in ['EF102815', 'EF100778']:
    r = run(['curl', '-s', '--max-time', '60',
             'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=%s&rettype=fasta&retmode=text' % acc], log=LOG)
    if r.returncode == 0 and r.stdout.startswith('>'):
        seq = ''.join(r.stdout.split('\n')[1:]).replace('\r', '').upper()
        ext16[acc] = seq
        LOG.write('fetched %s len=%d\n' % (acc, len(seq)))
    else:
        LOG.write('fetch %s failed\n' % acc)

# 16S 与各参考的比较：BLASTN（取最佳 HSP identity）
sixteen_rows = []
cand = [('smbu06_16S_copy1', s06_16s[0][1])]
cand += [('smbu08_16S_copy1', s08_16s[0][1])]
for ename, eseq in ext16.items():
    cand.append(('external_%s' % ename, eseq))
# 参考：194 16S 由 smbu06 16S 对 194 chr 的 BLAST 得到片段；直接比 smbu06/smbu08 16S 与 externals + 彼此
allseq = cand
fa_all = os.path.join(TMP, '16s_all.fa')
write_fasta(fa_all, [(n, s) for n, s in allseq])
db_all = os.path.join(TMP, '16s_db')
run([MAKEDB, '-in', fa_all, '-dbtype', 'nucl', '-out', db_all], log=None)
for qn, qs in allseq:
    r = run([BLASTN, '-query', fa_all, '-db', db_all, '-outfmt', '6 qseqid sseqid pident length qlen slen',
             '-evalue', '1e-10', '-max_target_seqs', '20'], log=None)
    hits = {}
    for line in r.stdout.strip().split('\n'):
        if not line:
            continue
        q, s, pid, alen, qlen, slen = line.split('\t')
        if q == qn and s != qn:
            hits[s] = (float(pid), int(alen), int(qlen), int(slen))
    for s, (pid, alen, qlen, slen) in hits.items():
        sixteen_rows.append((qn, s, pid, alen, qlen, slen))
# 去重：只保留双向最高的
with open(os.path.join(OUT, '16S_comparison.tsv'), 'w', encoding='utf-8') as f:
    f.write('query\tsubject\tpident\taln_len\tqlen\tslen\n')
    for r in sorted(set(sixteen_rows)):
        f.write('\t'.join(str(x) for x in r) + '\n')

# 16S vs 194 染色体：建立 194 染色体 DB 并定位 16S 拷贝
rfa194 = os.path.join(TMP, 'ref_194chr.fa')
write_fasta(rfa194, [('E_lactis_194_chr', refs['E_lactis_194']['chromosome'])])
db194 = os.path.join(TMP, 'db_194chr')
run([MAKEDB, '-in', rfa194, '-dbtype', 'nucl', '-out', db194], log=None, check=True)
r = run([BLASTN, '-query', fa_all, '-db', db194,
         '-outfmt', '6 qseqid pident length qstart qend sstart send slen', '-evalue', '1e-10',
         '-max_hsps', '20'], log=None)
with open(os.path.join(OUT, '16S_vs_194.tsv'), 'w', encoding='utf-8') as f:
    f.write('query\tpident\taln_len\tqstart\tqend\tsstart\tsend\t194_chr_len\n')
    for line in r.stdout.strip().split('\n'):
        if line:
            p = line.split('\t')
            if int(p[2]) >= 200:   # 只保留近全长命中
                f.write('\t'.join(p) + '\n')
LOG.write('16S vs 194 done\n')
LOG.write('DONE\n')
print('taxonomy_check done')
