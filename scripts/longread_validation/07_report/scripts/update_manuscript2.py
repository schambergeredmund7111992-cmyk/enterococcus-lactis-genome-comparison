# NOTE (release copy): absolute local paths were replaced with the placeholders
# <PROJECT_ROOT>, <DATA_ROOT>, <TOOLS_ROOT>, <HOME>, <PYTHON> etc. Adapt these at the
# top of the script (or via environment variables) before running. See README.md.
# -*- coding: utf-8 -*-
"""
update_manuscript2.py — 依据长读段验证结果更新论文（另存新文件，绝不覆盖原稿）
输入 : <WORKDIR>\\manuscript_fin_修订投稿版.docx
输出 : <WORKDIR>\\manuscript_fin_修订投稿版_长读段验证后.docx
       <WORKDIR>\\manuscript_fin_修订投稿版_长读段验证后.md
"""
import os
import shutil

import docx

SRC = r'<WORKDIR>\manuscript_fin_修订投稿版.docx'
DST = r'<WORKDIR>\manuscript_fin_修订投稿版_长读段验证后.docx'
MD = r'<WORKDIR>\manuscript_fin_修订投稿版_长读段验证后.md'

shutil.copyfile(SRC, DST)
doc = docx.Document(DST)


def replace_para(prefix, new_text, must=True):
    """按前缀定位段落并整体替换文本（保留段落样式）"""
    for p in doc.paragraphs:
        if p.text.strip().startswith(prefix):
            for r in list(p.runs):
                r.text = ''
            if p.runs:
                p.runs[0].text = new_text
            else:
                p.add_run(new_text)
            return True
    if must:
        raise SystemExit('paragraph not found: %s' % prefix[:50])
    return False


# ---- 图3 图注（原第38段） ----
replace_para(
    'Figure 3. smbu08 circular assembly component consistent',
    'Figure 3. smbu08 circular assembly component consistent with an excised tandem-dimer form of the '
    'phage-related element. (a) Structural model based on the att-flanked chromosomal form and the smbu08 '
    'junction. (b) Gene map of the 77.4-kb component. (c) Read-depth comparison. (d) Reverse-strand '
    'alignment to the smbu06 chromosome. (e) Long-read junction support (Supplementary Fig. S3; '
    'Supplementary Table S3): 1,883 Nanopore reads span the plasmid tandem boundary under pre-registered '
    'criteria, so the assembled component is not a junction-level assembly artefact; the tandem-dimer and '
    'monomer-circle interpretations are sequence-equivalent except for a 1-bp inter-unit indel that the '
    'available reads cannot resolve. Autonomous replication and induction are not established.')

# ---- §3.4 正文（原第40段） ----
replace_para(
    'In smbu08, the differential element was absent from the chromosome',
    'In smbu08, a 77,437-bp circular assembly component was present, while the delivered smbu08 chromosome '
    'assembly contained the excised attJ state. The component\'s first 38,719-bp unit was rotation-equivalent '
    'to the reverse complement of the smbu06 chromosomal differential region, and the second unit matched the '
    'first at 99.997% except for one 1-bp deletion. The component had 4.1–5.0-fold higher Illumina depth than '
    'the chromosome (1,777.6× versus 430.4×). Long-read validation (Supplementary Fig. S3; Supplementary '
    'Table S3) showed that 1,883 Nanopore reads (992 under strict criteria: single primary alignment, '
    'mapQ ≥ 20, ≥ 3-kb anchors on both sides) continuously span the plasmid tandem boundary, and that 41 '
    'reads place the chromosomal attJ state while rejecting the retained-state alternative '
    '(ΔAS 3.6–28.3 kb). The tandem-dimer and monomer-circle interpretations are sequence-equivalent except '
    'for the single 1-bp inter-unit indel, which the available reads cannot resolve. The component\'s high '
    'relative abundance (5.67-fold chromosome depth, bootstrap CI 5.51–5.83; tandem-duplication-corrected '
    'element abundance 4.0-fold, CI 3.94–4.05) does not establish autonomous replication.')

# ---- 新增段：混合群体（插到 §3.4 正文之后） ----
anchor = None
for p in doc.paragraphs:
    if p.text.strip().startswith('In smbu08, a 77,437-bp circular assembly component was present'):
        anchor = p
        break
if anchor is None:
    raise SystemExit('anchor paragraph not found')
new_p = anchor.insert_paragraph_before(
    'Long-read analysis additionally showed that the sequenced smbu08 population contains both chromosome '
    'states. Junction-spanning reads at the two state-specific boundaries were 257 (retained, chr06-like '
    'flank|attL) versus 46 (excised, chr08 flank|attJ), ≈5.6:1; an aligner-independent 20-mer adjacency test '
    'on the same reads gave 184 versus 35 mutually exclusive reads; and the same asymmetry was observed in '
    'independent Illumina and Nanopore depth profiles and validated against the pure-retained smbu06 control '
    '(195 versus 8 reads). The retained-state reads carry the smbu08-defining SNP allele (chromosome '
    'position 2,404,144 C: 462 C versus 5 T at that position), excluding contamination from smbu06/194. The '
    'delivered smbu08 chromosome assembly therefore represents the minority excised haplotype, whereas the '
    'population majority retains the element; both states co-exist with the high-copy excised-element '
    'plasmid.')
# 将新段移到 anchor 之后
anchor._p.addnext(new_p._p)

# ---- 讨论段（原第53段） ----
replace_para(
    'The differential region has a coherent phage-related architecture.',
    'The differential region has a coherent phage-related architecture. Its integrase, terminase, portal, '
    'head, tail, baseplate and lysis functions, together with the att-flanked chromosomal state, support '
    'classification as a PBSX-like phage-related element. In smbu08, the high-coverage 77.4-kb circular '
    'assembly component is consistent with a tandem-dimer form of the excised element, and long reads '
    'directly span its tandem boundary (1,883 reads; Supplementary Fig. S3). The same reads show that the '
    'sequenced population is a mixture in which the retained chromosome state predominates (~5.6:1), so '
    'excision is captured as an on-going, population-level process rather than a fixed genomic state; this '
    'is the principal genomic novelty of the comparison. However, the available sequence data do not '
    'establish induction, a replication mechanism, heritability or biological consequence. Those '
    'possibilities remain hypotheses.')

# ---- Limitations（原第57段） ----
replace_para(
    'Limitations. All conclusions are computational.',
    'Limitations. All conclusions are computational. Long reads support the attJ junction and span the '
    'plasmid tandem boundary, but they cannot distinguish a tandem dimer from a monomer circle at equal '
    'sequence abundance, and the 1-bp inter-unit difference was not independently resolved. The smbu08 '
    'population is a mixture of retained and excised chromosome states, and the delivered chromosome '
    'assembly represents the minority excised haplotype; population proportions are estimated from read '
    'counts and may be affected by coverage and length biases. The high-coverage circular component has '
    'not been experimentally induced or functionally characterised. Read mapping supplies sequence '
    'support but does not prove absolute sample purity. The public panel is availability-biased, so '
    'reported frequencies are descriptive rather than population estimates. Candidate bacteriocin loci '
    'have not been shown to be expressed, processed or active. Finally, static genome comparisons cannot '
    'infer the direction or timing of plasmid or cassette transfer.')

# ---- 补充材料说明（原第61段） ----
replace_para(
    'Supplementary Methods provide command-level detail.',
    'Supplementary Methods provide command-level detail. Supplementary Tables S1–S7 include assembly/read '
    'statistics, structural-variant annotation, long-read junction-support evidence (Table S3: junction '
    'counts, model competition, chromosome-state counts, depth ratios and copy-number estimates), '
    'candidate bacteriocin loci, public-panel metadata and distributions, and functional/safety context. '
    'Supplementary Figures S1–S4 contain read mapping, nisin controls, long-read junction support '
    '(Fig. S3, including the retained/excised population-state counts) and secondary annotation results.')

# ---- 摘要/结果段（原第27段）补一句 ----
replace_para(
    'The smbu06 assembly comprises a 2,676,906-bp chromosome',
    'The smbu06 assembly comprises a 2,676,906-bp chromosome and a 128,837-bp plasmid. The smbu08 assembly '
    'comprises a 2,638,187-bp chromosome (representing the minority excised haplotype of a mixed '
    'population; see Section 3.4), the same-length plasmid and a 77,437-bp circular assembly component. '
    'Reads mapped across 99.999–100% of every replicon at ≥1× coverage. Under the combined-reference '
    'mapping design, 100% of reads aligned and ≥99.99% were assigned to either project genome, providing '
    'no evidence of substantial non-target sequence contamination (Supplementary Fig. S1; Supplementary '
    'Table S1).')

doc.save(DST)
print('docx saved:', DST)

# ---- Markdown 版 ----
lines = []
for p in docx.Document(DST).paragraphs:
    t = p.text.rstrip()
    if not t:
        lines.append('')
        continue
    if t.startswith(('Figure ', 'Table ', 'Supplementary')):
        lines.append('**%s**' % t)
    elif len(t) < 90 and not t.endswith('.'):
        lines.append('## ' + t)
    else:
        lines.append(t)
open(MD, 'w', encoding='utf-8').write('# smbu06/smbu08 比较基因组学论文（长读段验证后修订版）\n\n' +
                                     '\n\n'.join(lines) + '\n')
print('markdown saved:', MD)
