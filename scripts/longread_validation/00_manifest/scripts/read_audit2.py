# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
read_audit2.py — 输入审计（本会话独立实现）：
  - Nanopore FASTQ/FASTA 定位、大小、read 数、总碱基、长度分布、N50、SHA-256
  - 组装/注释输入清单 + SHA-256
  - 工具版本
输出：00_manifest/nanopore_read_stats.tsv, input_inventory.tsv, tool_versions.txt
依赖：00_manifest/_read_lengths.tmp（gzip -dc | awk 预生成的全部 read 长度列）
"""
import hashlib
import os
import subprocess
import datetime

ATTJ = r'<PROJECT_ROOT>\longread_validation'
IN08 = r'<PROJECT_ROOT>\smbu08\smbu08'
IN06 = r'<PROJECT_ROOT>\smbu06\smbu06'
MAN = os.path.join(ATTJ, '00_manifest')

FQ = os.path.join(IN08, '1.Cleandata', 'smbu08.filtered_reads.fq.gz')
FA = os.path.join(IN08, '1.Cleandata', 'smbu08.filtered_reads.fa.gz')


def sha256_file(p, bs=1 << 22):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        while True:
            b = f.read(bs)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


L = [int(x) for x in open(os.path.join(MAN, '_read_lengths.tmp'))]
L.sort()
n = len(L)
tot = sum(L)
acc = 0
n50 = 0
for x in reversed(L):
    acc += x
    if acc >= tot / 2:
        n50 = x
        break


def q(p):
    return L[min(n - 1, int(n * p))]


lines = ['file\treads\tbases\tmean_len\tmedian_len\tmin_len\tmax_len\tread_n50\tp10\tp25\tp75\tp90\tp95\tp99\tlen_ge_7kb\tlen_ge_30kb']
lines.append('%s\t%d\t%d\t%.0f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d' % (
    os.path.basename(FQ), n, tot, tot / n, L[n // 2], L[0], L[-1], n50,
    q(.10), q(.25), q(.75), q(.90), q(.95), q(.99),
    sum(1 for x in L if x >= 7000), sum(1 for x in L if x >= 30000)))
lines.append('# length_source\tgzip -dc %s | awk NR%%4==2 -> %s' % (FQ, '_read_lengths.tmp'))

inv = ['item\tpath\tsize_bytes\tsha256']
inv.append(('nanopore_fq.gz\t%s\t%d\t%s' % (FQ, os.path.getsize(FQ), sha256_file(FQ))))
inv.append(('nanopore_fa.gz\t%s\t%d\t%s' % (FA, os.path.getsize(FA), sha256_file(FA))))
for nm, p in [
    ('smbu08_chromosome', os.path.join(IN08, '2.Assembly', 'smbu08.Chromosome.fasta')),
    ('smbu08_plasmid(Plasmid1+2)', os.path.join(IN08, '2.Assembly', 'smbu08.Plasmid.fasta')),
    ('smbu08_complete', os.path.join(IN08, '2.Assembly', 'smbu08.Complete.genome.fasta')),
    ('smbu08_genbank', os.path.join(IN08, '2.Assembly', 'smbu08.genome.gb')),
    ('smbu06_chromosome', os.path.join(IN06, '2.Assembly', 'smbu06.Chromosome.fasta')),
    ('smbu06_plasmid', os.path.join(IN06, '2.Assembly', 'smbu06.Plasmid.fasta')),
    ('smbu06_complete', os.path.join(IN06, '2.Assembly', 'smbu06.Complete.genome.fasta')),
]:
    inv.append('%s\t%s\t%d\t%s' % (nm, p, os.path.getsize(p), sha256_file(p)))

open(os.path.join(MAN, 'nanopore_read_stats.tsv'), 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
open(os.path.join(MAN, 'input_inventory.tsv'), 'w', encoding='utf-8').write('\n'.join(inv) + '\n')

mm2 = r'<TOOLS_ROOT>\minimap2-win\minimap2-2.31-r1302-windows-x86_64-ucrt64\minimap2.exe'
samtools = os.path.join(MAN, 'tools', 'bin', 'samtools.exe')
tv = ['attJ/串联二聚体长读段验证 — 工具与环境  %s' % datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')]
for nm, cmd in [('minimap2', [mm2, '--version']),
                ('blastn', [r'<TOOLS_ROOT>\ncbi-blast-2.17.0+\bin\blastn.exe', '-version']),
                ('python', [r'<PYTHON>', '--version'])]:
    r = subprocess.run(cmd, capture_output=True, text=True, encoding='utf-8', errors='replace')
    tv.append('%s: %s' % (nm, ((r.stdout or '') + (r.stderr or '')).strip().splitlines()[0]))
if os.path.exists(samtools):
    r = subprocess.run([samtools, '--version'], capture_output=True, text=True, encoding='utf-8', errors='replace')
    tv.append('samtools: %s (%s)' % (((r.stdout or '') + (r.stderr or '')).strip().splitlines()[0], samtools))
else:
    tv.append('samtools: 尚未就绪（MSYS2 mingw-w64-ucrt-x86_64-samtools 1.24-1 依赖闭包下载中）')
tv.append('说明：本机 Windows 无官方 samtools 分发；samtools 1.24-1 取自 MSYS2 官方镜像 mirror.msys2.org（ucrt64 包，含依赖闭包），'
          '位于 00_manifest/tools/bin。BAM 排序/索引/深度由该 samtools 完成。')
open(os.path.join(MAN, 'tool_versions.txt'), 'w', encoding='utf-8').write('\n'.join(tv) + '\n')
print('read_audit2 done')
for l in lines + inv[:3]:
    print(l)
