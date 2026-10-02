# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
p3lib.py — smbu06/smbu08 比较基因组项目公共库
约定：
  * 所有脚本通过 ASCII junction 路径访问文件（中文路径会让 BLAST+/MAFFT 静默失败；
    junction: <PROJECT_ROOT>\\{smbu06,smbu08,p3cmp,p3ws}）
  * 每次工具调用都记录命令行与返回码，写入调用方日志
"""
import hashlib
import os
import subprocess
import sys

# ---- ASCII 路径（junction）----
LINKS = r'<PROJECT_ROOT>'
WORK  = LINKS + r'\p3cmp'     # -> 论文3\smbu06_smbu08_比较基因组完整分析
IN06  = LINKS + r'\smbu06'    # -> 论文3\smbu06\smbu06
IN08  = LINKS + r'\smbu08'    # -> Downloads\smbu08\smbu08
P3WS  = LINKS + r'\p3ws'      # -> 论文3（整树）

# 真实（中文）路径，仅供文档与追溯
WORK_REAL = r'<PROJECT_ROOT>\comparative_analysis'

# ---- 工具 ----
TOOLS   = r'<TOOLS_ROOT>'
BLAST   = TOOLS + r'\ncbi-blast-2.17.0+\bin'
BLASTN  = BLAST + r'\blastn.exe'
TBLASTN = BLAST + r'\tblastn.exe'
BLASTP  = BLAST + r'\blastp.exe'
MAKEDB  = BLAST + r'\makeblastdb.exe'
MM2     = TOOLS + r'\minimap2-win\minimap2-2.31-r1302-windows-x86_64-ucrt64\minimap2.exe'
MAFFT   = TOOLS + r'\mafft-win\mafft.bat'
IQTREE  = TOOLS + r'\iqtree-2.3.6-Windows\bin\iqtree2.exe'
DATASETS = TOOLS + r'\datasets.exe'
PY      = r'<PYTHON>'

# 输入组装与注释（固定路径）
ASM06 = {
    'Chromosome1': IN06 + r'\smbu06\2.Assembly\smbu06.Chromosome.fasta',
    'Plasmid1':    IN06 + r'\smbu06\2.Assembly\smbu06.Plasmid.fasta',
    'complete':    IN06 + r'\smbu06\2.Assembly\smbu06.Complete.genome.fasta',
}
ASM08 = {
    'Chromosome1': IN08 + r'\smbu08\2.Assembly\smbu08.Chromosome.fasta',
    'Plasmid1':    IN08 + r'\smbu08\2.Assembly\smbu08.Plasmid.fasta',
    'complete':    IN08 + r'\smbu08\2.Assembly\smbu08.Complete.genome.fasta',
}
GFF06 = IN06 + r'\smbu06\3.Genome_Component\Gene_Predict\smbu06.Gene.gff'
GFF08 = IN08 + r'\smbu08\3.Genome_Component\Gene_Predict\smbu08.Gene.gff'
PEP06 = IN06 + r'\smbu06\3.Genome_Component\Gene_Predict\smbu06.Gene.pep.fasta'
PEP08 = IN08 + r'\smbu08\3.Genome_Component\Gene_Predict\smbu08.Gene.pep.fasta'
GBK06 = IN06 + r'\smbu06\2.Assembly\smbu06.genome.gb'
GBK08 = IN08 + r'\smbu08\2.Assembly\smbu08.genome.gb'


# ---- FASTA ----
def read_fasta(path):
    """返回 [(name, seq), ...] 保序；seq 大写。支持 .gz"""
    import gzip
    op = gzip.open if path.endswith('.gz') else open
    recs = []
    name, buf = None, []
    with op(path, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    recs.append((name, ''.join(buf).upper()))
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
    if name is not None:
        recs.append((name, ''.join(buf).upper()))
    return recs


def write_fasta(path, recs, width=70):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        for name, seq in recs:
            f.write('>' + name + '\n')
            for i in range(0, len(seq), width):
                f.write(seq[i:i + width] + '\n')


# ---- 环状序列规范化 ----
_COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def revcomp(s):
    return s.translate(_COMP)[::-1]


def min_rotation(s):
    """Booth 算法：字典序最小循环旋转 O(n)"""
    s = s + s
    n = len(s)
    f = [-1] * n
    k = 0
    for j in range(1, n):
        sj = s[j]
        i = f[j - k - 1]
        while i != -1 and sj != s[k + i + 1]:
            if sj < s[k + i + 1]:
                k = j - i - 1
            i = f[i]
        if sj != s[k + i + 1]:
            if sj < s[k]:
                k = j
            f[j - k] = -1
        else:
            f[j - k] = i + 1
    return s[k:k + n // 2]


def canonical_circular(seq):
    """双链环状序列的规范形式：两条链最小旋转中取字典序更小者"""
    a = min_rotation(seq)
    b = min_rotation(revcomp(seq))
    return min(a, b)


# ---- 哈希 ----
def sha256_file(path, blocksize=1 << 20):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        while True:
            b = f.read(blocksize)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def sha256_text(s):
    if isinstance(s, str):
        s = s.encode()
    return hashlib.sha256(s).hexdigest()


# ---- 运行与日志 ----
def run(cmd, log=None, check=False, cwd=None):
    """执行命令；cmd 为 list[str]。返回 CompletedProcess(text)。"""
    if log is not None:
        log.write('>>> ' + ' '.join(cmd) + '\n')
        log.flush()
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=cwd,
                       encoding='utf-8', errors='replace')
    if log is not None:
        log.write('<<< rc=%d\n' % r.returncode)
        if r.stdout:
            log.write(r.stdout[-4000:] + ('\n' if not r.stdout.endswith('\n') else ''))
        if r.stderr:
            log.write('[stderr] ' + r.stderr[-2000:] + '\n')
        log.flush()
    if check and r.returncode != 0:
        raise RuntimeError('command failed rc=%d: %s\n%s' % (r.returncode, ' '.join(cmd), r.stderr[-2000:]))
    return r


def open_log(path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    return open(path, 'a', encoding='utf-8')


def logpath(module, name):
    return os.path.join(WORK, module, 'logs', name)


def gc_content(s):
    s = s.upper()
    return (s.count('G') + s.count('C')) / len(s) * 100 if s else 0.0


# ---- 翻译（标准细菌遗传密码，TGA/TAA/TAG 终止） ----
_CODON = {}
_bases = 'TCAG'
_aas = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
_i = 0
for b1 in _bases:
    for b2 in _bases:
        for b3 in _bases:
            _CODON[b1 + b2 + b3] = _aas[_i]
            _i += 1


def translate_dna(s, to_stop=False):
    s = s.upper()
    n = len(s) - len(s) % 3
    aa = []
    for i in range(0, n, 3):
        c = _CODON.get(s[i:i + 3], 'X')
        if c == '*':
            if to_stop:
                break
            aa.append('*')
        else:
            aa.append(c)
    return ''.join(aa)
