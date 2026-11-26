import collections, statistics, os, csv
V = os.path.expanduser('~/sbria_anammox_shotgun/verificacao_itens_1_2/1e')
node = lambda s: '_'.join(s.split('_')[:2])
b9 = {node(l.split('\t')[0]) for l in open(V + '/bin009_copy.fasta.fai')}
ctrl = [(node(l.split('\t')[0]), int(l.split('\t')[1])) for l in open(V + '/ctrl_contigs.tsv')]
REP = {'NODE_17176', 'NODE_19120'}
L = collections.defaultdict(collections.Counter)
for l in open(V + '/partner_pairs.tsv'):
    c, m = l.rstrip('\n').split('\t'); L[node(c)][node(m)] += 1
def q(x):
    qs = statistics.quantiles(x, n=4, method='exclusive'); qi = statistics.quantiles(x, n=4, method='inclusive')
    return f'median {statistics.median(x):.3f}; Q1–Q3 {qs[0]:.3f}–{qs[2]:.3f} (exclusive) / {qi[0]:.3f}–{qi[2]:.3f} (inclusive); min {min(x):.3f}; max {max(x):.3f}'
out = open(V + '/control_22.tsv', 'w')
out.write('contig\tlength\tpairs_all\tto_bin009_all\tfrac_all\tpairs_to_repeats\tpairs_excl\tto_bin009_excl\tfrac_excl\n')
fa, fe = [], []
for c, ln in ctrl:
    d = L[c]; na = sum(d.values()); ba = sum(v for m, v in d.items() if m in b9)
    rp = sum(v for m, v in d.items() if m in REP); ne = na - rp; be = sum(v for m, v in d.items() if m in b9 and m not in REP)
    f1 = ba / na if na else float('nan'); f2 = be / ne if ne else float('nan')
    fa.append(f1); fe.append(f2)
    out.write(f'{c}\t{ln}\t{na}\t{ba}\t{f1:.3f}\t{rp}\t{ne}\t{be}\t{f2:.3f}\n')
out.close()
print(open(V + '/control_22.tsv').read())
print('repeats in bin.009?', REP & b9)
print('WITH repeat pairs   :', q(fa), f'| zeros {sum(f == 0 for f in fa)}')
print('WITHOUT repeat pairs:', q(fe), f'| zeros {sum(f == 0 for f in fe)}')
d = L['NODE_5312']; tot = sum(d.values()); tb = sum(v for m, v in d.items() if m in b9)
print(f'\nNODE_5312: {tot} inter-contig pairs; to bin.009 contigs: {tb}; to repeats: {sum(v for m,v in d.items() if m in REP)}')
for m, v in d.most_common(): print(f'   {m:12s} {v:4d}  {"bin.009" if m in b9 else "not bin.009"}')
four = ['NODE_573', 'NODE_5754', 'NODE_964', 'NODE_2986']
s4 = sum(d[m] for m in four); print(f'four named partners: {s4} of {tb} bin.009 pairs ({100*s4/tb:.1f}%); top-4 by count among bin.009 partners: {[m for m,_ in sorted(((m,v) for m,v in d.items() if m in b9), key=lambda x:-x[1])[:4]]}')
