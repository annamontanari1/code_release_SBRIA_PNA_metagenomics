import collections, os, csv
P = os.path.expanduser('~/sbria_anammox_shotgun'); V = P + '/verificacao_itens_1_2/1c'; T = V + '/taxdump'
cnt = collections.Counter(); st = collections.Counter()
with open(P + '/Kaiju/NGS1140_6S_trimmed_reads.kaiju') as f:
    for l in f:
        a = l.split('\t', 3); st[a[0]] += 1
        if a[0] == 'C': cnt[a[2].strip()] += 1
tot = sum(st.values()); C = st['C']
print(f'lines (read pairs): {tot:,}; C: {C:,} ({100*C/tot:.2f}%); U: {st["U"]:,}; distinct taxids: {len(cnt):,}')
parent, rank = {}, {}
for l in open(T + '/nodes.dmp'):
    x = l.split('\t|\t'); parent[x[0]] = x[1]; rank[x[0]] = x[2]
name = {}
for l in open(T + '/names.dmp'):
    x = l.split('\t|\t')
    if x[3].startswith('scientific name'): name[x[0]] = x[1]
merged = dict(l.split('\t|\t')[:2] for l in open(T + '/merged.dmp'))
merged = {k: v.strip('\t|\n') for k, v in merged.items()}
deleted = {l.split('\t')[0] for l in open(T + '/delnodes.dmp')}
def lineage(t):
    t = merged.get(t, t); out = []
    while t in parent and t not in out:
        out.append(t)
        if parent[t] == t: break
        t = parent[t]
    return out
targets = {}
for n in ['Candidatus Brocadia', 'Candidatus Kuenenia', 'Candidatus Brocadiaceae', 'Candidatus Brocadiales', 'Candidatus Brocadiia', 'Planctomycetota']:
    ids = [k for k, v in name.items() if v == n]; targets[n] = ids
print('target taxids (current NCBI):', {n: [(i, rank[i]) for i in ids] for n, ids in targets.items()})
tid = {n: set(ids) for n, ids in targets.items()}
sums = collections.Counter(); detail = collections.defaultdict(collections.Counter); missing = collections.Counter()
for t, c in cnt.items():
    lin = lineage(t)
    if not lin: missing[t] += c; continue
    for n, ids in tid.items():
        if ids & set(lin): sums[n] += c; detail[n][t] += c
print(f'taxids not found in current taxonomy: {len(missing)} ({sum(missing.values()):,} pairs; deleted: {sum(c for t,c in missing.items() if t in deleted):,})')
rows = []
for n in tid:
    s = sums[n]; rows.append((n, s, 100*s/C, 100*s/tot))
    print(f'{n:28s} {s:>9,} pairs  {100*s/C:7.4f}% of classified  {100*s/tot:7.4f}% of all filtered pairs')
    for t, c in detail[n].most_common(8): print(f'      {t:>8} {rank.get(merged.get(t,t),"?"):12s} {name.get(merged.get(t,t),"?"):45s} {c:>9,}')
with open(V + '/kaiju_brocadiaceae_counts.tsv', 'w') as f:
    f.write('taxon\tread_pairs\tpct_classified\tpct_all_filtered\n')
    for r in rows: f.write(f'{r[0]}\t{r[1]}\t{r[2]:.4f}\t{r[3]:.4f}\n')
# Brocadia taxids anywhere in output (any species/strain)
broc_desc = {k for k in parent if tid['Candidatus Brocadia'] & set(lineage(k))}
print(f'Brocadia descendant taxids in current NCBI: {len(broc_desc)}; present in .kaiju: {sum(1 for t in cnt if merged.get(t,t) in broc_desc)}')
# compare with previous tables
prev = {}
for r in csv.reader(open(P + '/results/review_checks/01_kaiju/kaiju_raw_counts.tsv'), delimiter='\t'):
    if r[0] == 'status': prev[r[1]] = int(r[2])
print('01_kaiju status counts match:', prev.get('C') == C and prev.get('U') == st['U'])
