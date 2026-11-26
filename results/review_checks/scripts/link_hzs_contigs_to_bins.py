#!/usr/bin/env python3
"""Liga os 10 contigs hzs/hdh aos bins por pares de leituras (Bowtie2 contra a montagem completa).

Entrada : results/contig_check/asm.sorted.bam (+ .bai), contigs->bins do DAS Tool (via base.py)
Saida   : results/review_checks/13_hzs_links/
  summary_links.tsv        1 linha por contig (leituras, pares inter-contig, ligacoes por bin)
  links_by_mate.tsv        cada contig vizinho (bin, tamanho, profundidade, ligacoes)
  control_bin009.tsv       controle positivo: contigs do proprio bin.009 (taxa de ligacao a bin.009)
  contig_context.tsv       GC e profundidade de cada contig vs bin.009
  SUMMARY.txt

Criterio: leitura primaria, mapeada, mate mapeado em OUTRO contig. MAPQ>=20 e a coluna conservadora
(leituras com MAPQ baixo vem de regioes repetidas e nao devem sustentar a ligacao).
"""
import os, sys, collections, statistics
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from base import *   # ROOT, RC, INFO, FINALBIN, node, bin_contigs
try:
    import pysam
except ImportError:
    sys.exit("pysam ausente. Instale: CONDA_SUBDIR=osx-arm64 conda install -y -n phylo --override-channels -c conda-forge -c bioconda pysam\n"
             "e rode com PATH incluindo ~/miniconda3/envs/phylo/bin (python do ambiente phylo).")

OUT = RC + '/13_hzs_links/'; os.makedirs(OUT, exist_ok=True)
BAM = ROOT + '/results/contig_check/asm.sorted.bam'
TARGET_BIN = 'bin.009'
CONTIGS = ['NODE_5668', 'NODE_5312', 'NODE_7170', 'NODE_17695', 'NODE_29475',
           'NODE_24575', 'NODE_20914', 'NODE_23806', 'NODE_47442', 'NODE_24842']
GENES = {'NODE_5668': 'hzsABC', 'NODE_5312': 'hzsABC', 'NODE_7170': 'hzsABC', 'NODE_17695': 'hzsAB(C)',
         'NODE_29475': 'hzsA', 'NODE_24575': 'hzsB', 'NODE_20914': 'hdh', 'NODE_23806': 'hdh',
         'NODE_47442': 'hdh', 'NODE_24842': 'hdh'}
MAPQ_MIN = 20

bam = pysam.AlignmentFile(BAM)
n2 = {r: node(r) for r in bam.references}

def bin_of(c): return FINALBIN.get(c, 'unbinned')

def scan(c):
    """conta leituras de c e ligacoes a outros contigs."""
    allL, q20L = collections.Counter(), collections.Counter()
    n = unm = same = 0
    for r in bam.fetch(INFO[c]['full']):
        if r.is_secondary or r.is_supplementary or r.is_unmapped: continue
        n += 1
        if r.mate_is_unmapped: unm += 1; continue
        m = n2[r.next_reference_name]
        if m == c: same += 1; continue
        allL[m] += 1
        if r.mapping_quality >= MAPQ_MIN: q20L[m] += 1
    return n, unm, same, allL, q20L

def by_bin(cnt):
    d = collections.Counter()
    for m, v in cnt.items(): d[bin_of(m)] += v
    return d

# ---- alvos
rows, mates = [], []
for c in CONTIGS:
    n, unm, same, allL, q20L = scan(c)
    b_all, b_q20 = by_bin(allL), by_bin(q20L)
    tot_q20 = sum(q20L.values())
    rows.append(dict(contig=c, genes=GENES[c], length=INFO[c]['len'], reads=n, mate_unmapped=unm, mate_same=same,
                     inter_all=sum(allL.values()), inter_q20=tot_q20,
                     to_bin009_all=b_all.get(TARGET_BIN, 0), to_bin009_q20=b_q20.get(TARGET_BIN, 0),
                     to_unbinned_q20=b_q20.get('unbinned', 0),
                     to_other_bins_q20='; '.join(f'{b}:{v}' for b, v in b_q20.most_common() if b not in (TARGET_BIN, 'unbinned')) or '-'))
    for m, v in allL.most_common():
        mates.append((c, m, bin_of(m), INFO[m]['len'], round(INFO[m]['depth'], 2), v, q20L.get(m, 0)))

# ---- controle positivo: contigs do bin.009 (1-5 kb, como os alvos) -> fracao de ligacoes q20 que vao ao bin.009
ctrl = []
for c in sorted(bin_contigs(TARGET_BIN), key=lambda x: INFO[x]['len']):
    if not (1000 <= INFO[c]['len'] <= 6000): continue
    n, unm, same, allL, q20L = scan(c)
    b = by_bin(q20L); tot = sum(q20L.values())
    ctrl.append((c, INFO[c]['len'], n, tot, b.get(TARGET_BIN, 0), round(b.get(TARGET_BIN, 0) / tot, 3) if tot else 'NA'))
ctrl_fracs = [x[5] for x in ctrl if x[5] != 'NA']
ctrl_links = [x[4] for x in ctrl]

# ---- contexto GC / profundidade
b9 = bin_contigs(TARGET_BIN)
gc9 = [INFO[c]['gc'] for c in b9]; dp9 = [INFO[c]['depth'] for c in b9]
gc_m, gc_s = statistics.mean(gc9), statistics.stdev(gc9); dp_med = statistics.median(dp9)
ctx = [(c, INFO[c]['len'], INFO[c]['gc'], round((INFO[c]['gc'] - gc_m) / gc_s, 2), round(INFO[c]['depth'], 2),
        round(INFO[c]['depth'] / dp_med, 2), bin_of(c)) for c in CONTIGS]

# ---- escrita
def w(path, header, data):
    with open(OUT + path, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for r in data: f.write('\t'.join(map(str, r)) + '\n')
H = ['contig', 'genes', 'length', 'reads', 'mate_unmapped', 'mate_same_contig', 'inter_contig_all', 'inter_contig_q20',
     'to_bin009_all', 'to_bin009_q20', 'to_unbinned_q20', 'to_other_bins_q20']
w('summary_links.tsv', H, [[r[k] for k in ['contig', 'genes', 'length', 'reads', 'mate_unmapped', 'mate_same', 'inter_all', 'inter_q20',
                                           'to_bin009_all', 'to_bin009_q20', 'to_unbinned_q20', 'to_other_bins_q20']] for r in rows])
w('links_by_mate.tsv', ['contig', 'mate_contig', 'mate_bin', 'mate_len', 'mate_depth', 'links_all', 'links_q20'], mates)
w('control_bin009.tsv', ['bin009_contig', 'length', 'reads', 'inter_q20', 'q20_to_bin009', 'frac_to_bin009'], ctrl)
w('contig_context.tsv', ['contig', 'length', 'GC_pct', 'GC_z_vs_bin009', 'depth_x', 'depth_ratio_vs_bin009_median', 'final_bin'], ctx)

with open(OUT + 'SUMMARY.txt', 'w') as f:
    f.write(f'bin.009: {len(b9)} contigs; GC {gc_m:.2f} +- {gc_s:.2f}%; profundidade mediana {dp_med:.2f}x\n')
    f.write(f'Controle (contigs do bin.009 de 1-6 kb, n={len(ctrl)}): fracao mediana de ligacoes q20 que vao ao bin.009 = '
            f'{statistics.median(ctrl_fracs):.2f}; ligacoes q20 ao bin.009 por contig: mediana {statistics.median(ctrl_links):.0f} '
            f'(min {min(ctrl_links)}, max {max(ctrl_links)})\n\n')
    for r in rows:
        f.write(f"{r['contig']} ({r['genes']}, {r['length']} bp, {r['reads']} leituras): pares inter-contig {r['inter_all']} (q20: {r['inter_q20']}); "
                f"-> bin.009: {r['to_bin009_all']} (q20: {r['to_bin009_q20']}); -> sem bin (q20): {r['to_unbinned_q20']}; outros bins (q20): {r['to_other_bins_q20']}\n")
    f.write('\nVizinhos q20 que mais ligam (por contig):\n')
    for c in CONTIGS:
        top = [m for m in mates if m[0] == c and m[6] > 0][:4]
        f.write(f'  {c}: ' + (', '.join(f'{m[1]} [{m[2]}, {m[5]}/{m[6]}]' for m in top) or 'nenhum com MAPQ>=20') + '\n')
print(open(OUT + 'SUMMARY.txt').read())
print('Saidas em', OUT)
