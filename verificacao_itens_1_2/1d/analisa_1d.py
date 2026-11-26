import csv, os, collections
P = os.path.expanduser('~/sbria_anammox_shotgun'); V = P + '/verificacao_itens_1_2/1d'
s2b = {}
for l in open(P + '/das_tool_result (1)/das_tool_output_dir_DASTool_scaffolds2bin.txt'):
    c, b = l.rstrip('\n').split('\t')[:2]; s2b[c] = b.replace('.fasta', '').replace('.fa', '')
bins = sorted(set(s2b.values()))
cov = {r['#rname']: int(r['numreads']) for r in csv.DictReader(open(V + '/coverage.tsv'), delimiter='\t')}
prim, q20 = {}, {}
for l in open(V + '/primary_per_contig.tsv'):
    c, n, q = l.split('\t'); prim[c] = int(n); q20[c] = int(q)
b9c = {l[1:].split()[0] for l in open(P + '/das_tool_result (1)/das_tool_output_dir_DASTool_bins/bin.009.fasta') if l.startswith('>')}
src = {s2b[c] for c in b9c}; print('bin.009.fasta contigs:', len(b9c), '| source bin in scaffolds2bin:', src, '| that bin has', sum(1 for b in s2b.values() if b in src), 'contigs')
for c in b9c: s2b[c] = 'bin.009'
assert set(s2b) <= set(cov), 'contig names in scaffolds2bin not found in BAM'
TOTAL = 44906144; ALN = sum(cov.values())
print(f'bins: {len(bins)}; binned contigs: {len(s2b)}; contigs in BAM: {len(cov)}')
print(f'filtered reads: {TOTAL:,} (= 2 x {TOTAL//2:,} pairs); aligned (numreads sum): {ALN:,} = {100*ALN/TOTAL:.2f}%')
print(f'aligned (primary per-contig sum): {sum(prim.values()):,}; MAPQ>=20: {sum(q20.values()):,} ({100*sum(q20.values())/TOTAL:.2f}% of filtered)')
def tot(d, which): return sum(d[c] for c, b in s2b.items() if which is None or b == which)
rows = []
for label, d, aln in [('published (samtools coverage numreads)', cov, ALN), ('primary mapped (-F 0x904)', prim, sum(prim.values())), ('primary, MAPQ>=20', q20, sum(q20.values()))]:
    a25 = tot(d, None); b9 = tot(d, 'bin.009')
    rows.append((label, a25, b9, aln))
    print(f'\n[{label}]\n  25 bins: {a25:,} reads = {100*a25/aln:.2f}% of aligned = {100*a25/TOTAL:.2f}% of filtered')
    print(f'  bin.009: {b9:,} reads = {100*b9/TOTAL:.2f}% of filtered = {100*b9/aln:.2f}% of aligned = {100*b9/a25:.2f}% of the 25-bin reads')
print('\nbin.009 contigs:', sum(1 for b in s2b.values() if b == 'bin.009'))
with open(V + '/bin009_recruitment.tsv', 'w') as f:
    f.write('method\treads_25bins\treads_bin009\taligned_total\tpct_bin009_of_filtered\tpct_bin009_of_aligned\tpct_bin009_of_25bins\n')
    for l, a25, b9, aln in rows: f.write(f'{l}\t{a25}\t{b9}\t{aln}\t{100*b9/TOTAL:.4f}\t{100*b9/aln:.4f}\t{100*b9/a25:.4f}\n')
