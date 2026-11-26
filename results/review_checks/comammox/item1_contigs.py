"""Item 1: NODE_2578 e NODE_13511 — métricas, TNF vs bin.007/bin.004, binners, ordem gênica, pares de reads."""
import sys, os, random, itertools, collections
sys.path.insert(0, '/private/tmp/claude-501/-Users-annamontanari/db506465-ed69-4ea1-9e3d-4df8c88155e2/scratchpad/scripts')
from base import *
import numpy as np, pysam
OUT = RC + '/comammox/'
C = ['NODE_2578', 'NODE_13511']
need = set(C) | set(bin_contigs('bin.007')) | set(bin_contigs('bin.004'))
SEQ = {}; k = None
for l in open(ROOT + '/results/targeted_search/bin_assemblies/NGS1140_6S_assembly.fa'):
    if l[0] == '>': n = node(l[1:]); k = n if n in need else None; SEQ.setdefault(k, []) if k else None
    elif k: SEQ[k].append(l.strip())
SEQ = {a: ''.join(b).upper() for a, b in SEQ.items()}
comp = str.maketrans('ACGT', 'TGCA'); km = [''.join(p) for p in itertools.product('ACGT', repeat=4)]
can = sorted(set(min(x, x.translate(comp)[::-1]) for x in km)); IDX = {x: can.index(min(x, x.translate(comp)[::-1])) for x in km}
def tnf(s):
    v = np.zeros(len(can))
    for i in range(len(s) - 3):
        j = IDX.get(s[i:i + 4])
        if j is not None: v[j] += 1
    return v / v.sum()
random.seed(1)
cent, pools = {}, {}
for b in ['bin.007', 'bin.004']:
    cs = [c for c in bin_contigs(b) if c not in C]
    v = sum(tnf(SEQ[c]) * len(SEQ[c]) for c in cs); cent[b] = v / v.sum(); pools[b] = [SEQ[c] for c in cs]
def null(b, L, n=500):
    xs = pools[b]; w = [len(x) for x in xs]; d = []
    while len(d) < n:
        x = random.choices(xs, weights=w)[0]
        if len(x) >= L: i = random.randrange(len(x) - L + 1); d.append(np.linalg.norm(tnf(x[i:i + L]) - cent[b]))
    return np.array(d)
gcs = {b: np.array([INFO[c]['gc'] for c in bin_contigs(b)]) for b in cent}
dep = {b: np.array([INFO[c]['depth'] for c in bin_contigs(b)]) for b in cent}
with open(OUT + 'item1_contig_metrics_TNF.tsv', 'w') as w:
    w.write('contig\tlength_bp\tdepth_x\tGC_pct\tfinal_bin\tDAS_Tool_retained\tbinners\t'
            'TNF_dist_bin007\tTNF_null_P95_bin007\tTNF_pctile_bin007\tTNF_dist_bin004\tTNF_null_P95_bin004\tTNF_pctile_bin004\t'
            'bin007_GC_mean±sd\tbin004_GC_mean±sd\tbin007_depth_median\tbin004_depth_median\n')
    for c in C:
        I = INFO[c]; row = [c, I['len'], round(I['depth'], 2), I['gc'], FINALBIN.get(c, 'unbinned'), 'yes' if c in FINAL else 'no', binners(c)]
        for b in ['bin.007', 'bin.004']:
            d = np.linalg.norm(tnf(SEQ[c]) - cent[b]); nd = null(b, len(SEQ[c]))
            row += [round(d, 4), round(np.percentile(nd, 95), 4), round(100 * (nd < d).mean(), 1)]
        row += [f"{gcs['bin.007'].mean():.2f}±{gcs['bin.007'].std(ddof=1):.2f}", f"{gcs['bin.004'].mean():.2f}±{gcs['bin.004'].std(ddof=1):.2f}",
                round(float(np.median(dep['bin.007'])), 2), round(float(np.median(dep['bin.004'])), 2)]
        w.write('\t'.join(map(str, row)) + '\n')
# pares de reads
bam = pysam.AlignmentFile(ROOT + '/results/contig_check/asm.sorted.bam')
n2 = {r: node(r) for r in bam.references}
with open(OUT + 'item1_readpair_links.tsv', 'w') as w, open(OUT + 'item1_readpair_links_by_bin.tsv', 'w') as wb:
    w.write('contig\tmate_contig\tmate_final_bin\tmate_len\tmate_depth\tlinks_all\tlinks_MAPQ>=20\n')
    wb.write('contig\tprimary_reads\tmate_unmapped\tmate_same_contig\tmate_other_contig\tlinks_by_final_bin(all MAPQ; MAPQ>=20)\n')
    for c in C:
        L = collections.Counter(); L20 = collections.Counter(); n = un = same = 0
        for r in bam.fetch(INFO[c]['full']):
            if r.is_secondary or r.is_supplementary or r.is_unmapped: continue
            n += 1
            if r.mate_is_unmapped: un += 1; continue
            m = n2[r.next_reference_name]
            if m == c: same += 1; continue
            L[m] += 1
            if r.mapping_quality >= 20: L20[m] += 1
        for m, v in L.most_common(): w.write(f"{c}\t{m}\t{FINALBIN.get(m, 'unbinned')}\t{INFO[m]['len']}\t{INFO[m]['depth']:.2f}\t{v}\t{L20[m]}\n")
        byb = collections.Counter(); byb20 = collections.Counter()
        for m, v in L.items(): byb[FINALBIN.get(m, 'unbinned')] += v; byb20[FINALBIN.get(m, 'unbinned')] += L20[m]
        wb.write(f"{c}\t{n}\t{un}\t{same}\t{sum(L.values())}\t" + '; '.join(f'{b}: {v} ({byb20[b]})' for b, v in byb.most_common()) + '\n')
