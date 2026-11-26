#!/usr/bin/env python3
"""Secao 18 do README (23/10/2026): pares de leitura, profundidade e GC dos contigs candidatos (base da figura composta).

Regra dos pares (= link_hzs_contigs_to_bins.py e registros_15.py): leitura primaria e mapeada (sem secondary/supplementary),
mate mapeado em OUTRO contig; coluna principal = MAPQ >= 20 da leitura. BAM: results/contig_check/asm.sorted.bam (Bowtie2 2.5.5).
Profundidade: meandepth de results/contig_check/asm_coverage.tsv (samtools coverage). GC: 00_base/contig_len_gc.tsv.
Bins: DAS Tool scaffolds2bin.txt (via base.py).

Uso: python attribution_18.py task1   -> 18_attribution/task1_reproduction.tsv (valores publicados x recalculados)
     python attribution_18.py task2   -> 18_attribution/candidates_links.tsv, candidates_links_by_bin.tsv (todos os candidatos)
     python attribution_18.py panelA  -> 18_attribution/panelA_contigs.tsv (contigs >= 5 kb + todos os de bin.009/bin.007 + candidatos)
"""
import os, sys, collections, statistics
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from base import *   # ROOT, RC, INFO (len, gc, depth), FINALBIN, node
import pysam

OUT = RC + '/18_attribution/'; os.makedirs(OUT, exist_ok=True)
MAPQ_MIN = 20
HZS = ['NODE_5312', 'NODE_23806', 'NODE_5668', 'NODE_20914', 'NODE_7170', 'NODE_17695', 'NODE_29475', 'NODE_24575', 'NODE_47442', 'NODE_24842']
OTHER = ['NODE_13511', 'NODE_2578', 'NODE_4298', 'NODE_1742', 'NODE_696', 'NODE_2432', 'NODE_4969']
TARGET = {**{c: 'bin.009' for c in HZS + ['NODE_4298', 'NODE_1742', 'NODE_696']}, 'NODE_13511': 'bin.007', 'NODE_2578': 'bin.007',
          'NODE_2432': 'bin.004', 'NODE_4969': '-'}
bam = pysam.AlignmentFile(ROOT + '/results/contig_check/asm.sorted.bam')
n2c = {r: node(r) for r in bam.references}

def links(c):
    tot = unm = same = 0; L, L20 = collections.Counter(), collections.Counter()
    for r in bam.fetch(INFO[c]['full']):
        if r.is_secondary or r.is_supplementary or r.is_unmapped: continue
        tot += 1
        if r.mate_is_unmapped: unm += 1; continue
        m = n2c[r.next_reference_name]
        if m == c: same += 1; continue
        L[m] += 1
        if r.mapping_quality >= MAPQ_MIN: L20[m] += 1
    by, by20 = collections.Counter(), collections.Counter()
    for m in L: by[FINALBIN.get(m, 'unbinned')] += L[m]; by20[FINALBIN.get(m, 'unbinned')] += L20[m]
    return dict(reads=tot, mate_unmapped=unm, mate_same=same, inter_all=sum(L.values()), inter_q20=sum(L20.values()), L=L, L20=L20, by=by, by20=by20)

b9 = [c for c, b in FINALBIN.items() if b == 'bin.009']
gc9 = [INFO[c]['gc'] for c in b9]; d9 = [INFO[c]['depth'] for c in b9]
def wmedian(vals, w):
    o = sorted(zip(vals, w)); half = sum(w) / 2; acc = 0
    for v, x in o:
        acc += x
        if acc >= half: return v

if sys.argv[1:] == ['task1']:
    rows = []
    def add(item, pub, new): rows.append((item, pub, new, 'OK' if str(pub) == str(new) else 'DIVERGE')); print(f'{item}: publicado {pub} | recalculado {new} | {rows[-1][3]}')
    for c, pub in [('NODE_5312', '151/160'), ('NODE_23806', '151/153'), ('NODE_5668', '5/11')]:
        x = links(c); add(f'{c} pares q20 bin.009/inter-contig', pub, f"{x['by20']['bin.009']}/{x['inter_q20']}")
    x = links('NODE_20914'); add('NODE_20914 pares q20 com bin.009', 3, x['by20']['bin.009'])
    print(f"   (NODE_20914 detalhe q20: inter-contig {x['inter_q20']}; por bin {dict(x['by20'])}; por contig do bin.009 "
          f"{ {m: v for m, v in x['L20'].items() if FINALBIN.get(m) == 'bin.009'} })")
    x = links('NODE_13511')
    add('NODE_13511 pares q20 com bin.007', 46, x['by20']['bin.007']); add('NODE_13511 pares q20 com nao binados', 14, x['by20']['unbinned'])
    add('NODE_13511 pares q20 com outros bins', 0, sum(v for b, v in x['by20'].items() if b not in ('bin.007', 'unbinned')))
    add('NODE_5312 profundidade (x)', 54.6, round(INFO['NODE_5312']['depth'], 1)); add('NODE_23806 profundidade (x)', 58.6, round(INFO['NODE_23806']['depth'], 1))
    add('bin.009 profundidade, mediana ponderada por comprimento (x)', 26.7, round(wmedian(d9, [INFO[c]['len'] for c in b9]), 1))
    add('NODE_5312 GC (%)', 51.2, round(INFO['NODE_5312']['gc'], 1)); add('NODE_23806 GC (%)', 48.6, round(INFO['NODE_23806']['gc'], 1))
    add('bin.009 GC media +- DP (%, nao ponderada)', '42.5 +- 1.9', f'{statistics.mean(gc9):.1f} +- {statistics.stdev(gc9):.1f}')
    with open(OUT + 'task1_reproduction.tsv', 'w') as f:
        f.write('item\tpublished\trecomputed\tstatus\n'); [f.write('\t'.join(map(str, r)) + '\n') for r in rows]
    print('divergencias:', sum(r[3] == 'DIVERGE' for r in rows))

if sys.argv[1:] == ['task2']:
    with open(OUT + 'candidates_links.tsv', 'w') as f, open(OUT + 'candidates_links_by_bin.tsv', 'w') as g:
        f.write('contig\tgroup\ttarget_bin\tfinal_bin\tlength_bp\tdepth_x\tGC_pct\treads\tmate_unmapped\tmate_same\tinter_all\tinter_q20\t'
                'q20_to_target\tq20_to_other_bins\tq20_to_unbinned\tfrac_target\tfrac_other\tfrac_unbinned\n')
        g.write('contig\tbin\tlinks_all\tlinks_q20\ttop_q20_contigs\n')
        for c in HZS + OTHER:
            x = links(c); t = TARGET[c]; q = x['inter_q20']
            tgt = x['by20'][t] if t != '-' else 0; un = x['by20']['unbinned']; oth = q - tgt - un
            grp = 'hzs/hdh' if c in HZS else {'NODE_13511': 'amoAB', 'NODE_2578': 'amoC', 'NODE_4298': 'nxrA', 'NODE_1742': 'nxrA',
                                              'NODE_696': 'nirS', 'NODE_2432': 'nxrA', 'NODE_4969': 'amoCAB'}[c]
            fr = lambda v: f'{v / q:.3f}' if q else 'NA'
            f.write(f"{c}\t{grp}\t{t}\t{FINALBIN.get(c, 'unbinned')}\t{INFO[c]['len']}\t{INFO[c]['depth']:.1f}\t{INFO[c]['gc']:.1f}\t{x['reads']}\t"
                    f"{x['mate_unmapped']}\t{x['mate_same']}\t{x['inter_all']}\t{q}\t{tgt}\t{oth}\t{un}\t{fr(tgt)}\t{fr(oth)}\t{fr(un)}\n")
            for b in sorted(x['by'], key=lambda b: -x['by20'][b]):
                top = ', '.join(f'{m}:{v}' for m, v in x['L20'].most_common() if FINALBIN.get(m, 'unbinned') == b)
                g.write(f"{c}\t{b}\t{x['by'][b]}\t{x['by20'][b]}\t{top}\n")
    print(open(OUT + 'candidates_links.tsv').read())

if sys.argv[1:] == ['panelA']:
    with open(OUT + 'panelA_contigs.tsv', 'w') as f:
        f.write('contig\tlength_bp\tGC_pct\tdepth_x\tfinal_bin\tclass\n')
        for c, d in INFO.items():
            b = FINALBIN.get(c, 'unbinned')
            cl = 'candidate' if c in TARGET else (b if b in ('bin.009', 'bin.007') else ('background' if d['len'] >= 5000 else None))
            if cl: f.write(f"{c}\t{d['len']}\t{d['gc']:.2f}\t{d['depth']:.3f}\t{b}\t{cl}\n")
    import collections as _c
    print(_c.Counter(l.rstrip().split('\t')[-1] for l in open(OUT + 'panelA_contigs.tsv') if not l.startswith('contig')))

# ---- v2 (23/10/2026, seção 19): NODE_16231 e pares com contigs repetitivos (NODE_17176, 148×; NODE_19120, 472×) fora da contagem
REPEATS = {'NODE_17176', 'NODE_19120'}
if sys.argv[1:] == ['v2']:
    TARGET['NODE_16231'] = 'bin.009'
    with open(OUT + 'candidates_links_v2.tsv', 'w') as f, open(OUT + 'NODE_16231_partners.tsv', 'w') as g:
        f.write('contig\ttarget_bin\tfinal_bin\tlength_bp\tdepth_x\tGC_pct\tinter_q20_all\tq20_to_repeats\tinter_q20\tq20_to_target\t'
                'q20_to_other_bins\tq20_to_unbinned\n')
        g.write('partner\tpairs_q20\tpairs_all\tlength_bp\tdepth_x\tGC_pct\tfinal_bin\trepeat_excluded\n')
        for c in HZS + ['NODE_13511', 'NODE_16231', 'NODE_4298', 'NODE_1742', 'NODE_696']:
            x = links(c); t = TARGET[c]
            rep = sum(v for m, v in x['L20'].items() if m in REPEATS)
            L = {m: v for m, v in x['L20'].items() if m not in REPEATS}
            q = sum(L.values()); tgt = sum(v for m, v in L.items() if FINALBIN.get(m) == t)
            un = sum(v for m, v in L.items() if m not in FINALBIN); oth = q - tgt - un
            f.write(f"{c}\t{t}\t{FINALBIN.get(c, 'unbinned')}\t{INFO[c]['len']}\t{INFO[c]['depth']:.1f}\t{INFO[c]['gc']:.1f}\t"
                    f"{x['inter_q20']}\t{rep}\t{q}\t{tgt}\t{oth}\t{un}\n")
            if c == 'NODE_16231':
                for m, v in x['L20'].most_common():
                    g.write(f"{m}\t{v}\t{x['L'][m]}\t{INFO[m]['len']}\t{INFO[m]['depth']:.1f}\t{INFO[m]['gc']:.1f}\t"
                            f"{FINALBIN.get(m, 'unbinned')}\t{'yes' if m in REPEATS else 'no'}\n")
    print(open(OUT + 'candidates_links_v2.tsv').read()); print(open(OUT + 'NODE_16231_partners.tsv').read())

if sys.argv[1:] == ['panelA_v2']:   # como panelA, com NODE_16231 entre os candidatos (v2 da figura)
    TARGET['NODE_16231'] = 'bin.009'
    with open(OUT + 'panelA_contigs_v2.tsv', 'w') as f:
        f.write('contig\tlength_bp\tGC_pct\tdepth_x\tfinal_bin\tclass\n')
        for c, d in INFO.items():
            b = FINALBIN.get(c, 'unbinned')
            cl = 'candidate' if c in TARGET else (b if b in ('bin.009', 'bin.007') else ('background' if d['len'] >= 5000 else None))
            if cl: f.write(f"{c}\t{d['len']}\t{d['gc']:.2f}\t{d['depth']:.3f}\t{b}\t{cl}\n")
