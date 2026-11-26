import csv, hashlib, os, re, collections
V = os.path.expanduser('~/sbria_anammox_shotgun/verificacao_itens_1_2/1a')
OLD = os.path.expanduser('~/sbria_anammox_shotgun/results/review_checks/35_skani_bin009')
acc_of = lambda p: re.search(r'(GC[AF]_\d+\.\d+)', os.path.basename(p)).group(1)
refs = [l.strip() for l in open(f'{V}/ref_list.txt')]
old = {r['Assembly Accession']: r for r in csv.DictReader(open(f'{OLD}/brocadia_genomes.tsv'), delimiter='\t')}
new = {r['Assembly Accession']: r for r in csv.DictReader(open(f'{V}/brocadia_records_20261123.tsv'), delimiter='\t')}
print(f'records 24/10 (brocadia_genomes.tsv): {len(old)}; NCBI today: {len(new)}; files in genomes_raw: {len(refs)}')
print('files == 24/10 list:', sorted(acc_of(p) for p in refs) == sorted(old))
print('new since 24/10:', sorted(set(new) - set(old)), '| removed:', sorted(set(old) - set(new)))
# (ii-a) pairing by NCBI paired accession (today's metadata, restricted to the 67)
def key_pair(a):
    p = new.get(a, {}).get('Assembly Paired Assembly Accession') or old[a].get('Assembly Paired Assembly Accession') or ''
    ids = [a] + ([p] if p else [])
    return min(re.sub(r'^GC[AF]_', '', i) for i in ids)
pairs = collections.defaultdict(list)
for a in old: pairs[key_pair(a)].append(a)
dup_pairs = {k: v for k, v in pairs.items() if len(v) > 1}
print(f'unique by paired accession: {len(pairs)}; pairs with 2 records: {sorted(dup_pairs.values())}')
# (ii-b) content checksum
def gmd5(p):
    seqs, cur = [], []
    for l in open(p):
        if l.startswith('>'):
            if cur: seqs.append(hashlib.md5(''.join(cur).upper().encode()).hexdigest()); cur = []
        else: cur.append(l.strip())
    if cur: seqs.append(hashlib.md5(''.join(cur).upper().encode()).hexdigest())
    return hashlib.md5(''.join(sorted(seqs)).encode()).hexdigest()
chk = {acc_of(p): gmd5(p) for p in refs}
byc = collections.defaultdict(list)
for a, c in chk.items(): byc[c].append(a)
print(f'unique by sequence checksum: {len(byc)}; identical groups: {sorted(sorted(v) for v in byc.values() if len(v) > 1)}')
uid = {a: chk[a] for a in chk}   # unique-assembly key = content checksum
# (iii/iv) skani rerun
sk = list(csv.DictReader(open(f'{V}/skani_rerun.tsv'), delimiter='\t'))
for r in sk: r['acc'] = acc_of(r['Ref_file'])
rep = {r['acc'] for r in sk}
print(f'skani reported: {len(rep)} records = {len({uid[a] for a in rep})} unique assemblies')
miss = sorted(set(old) - rep)
print(f'not reported: {len(miss)} records = {len({uid[a] for a in miss} - {uid[a] for a in rep})} unique (absent entirely)')
for a in miss: print('   ', a, old[a]['Organism Name'], old[a]['Organism Infraspecific Names Strain'] or old[a]['Organism Infraspecific Names Isolate'], old[a]['Assembly Level'])
# compare with original run
o = {acc_of(r['Ref_file']): r for r in csv.DictReader(open(f'{OLD}/skani_bin009_vs_brocadia.tsv'), delimiter='\t')}
same = all(o[a]['ANI'] == r['ANI'] and o[a]['Align_fraction_ref'] == r['Align_fraction_ref'] and o[a]['Align_fraction_query'] == r['Align_fraction_query'] for a, r in ((r['acc'], r) for r in sk)) and set(o) == rep
print('rerun identical to 24/10 run (ANI and both AF, same set):', same)
# (v) top 20
sk.sort(key=lambda r: -float(r['ANI']))
with open(f'{V}/top20.tsv', 'w') as f:
    f.write('rank\taccession\torganism\tstrain_isolate\tlevel\trelease_date\tANI\tAF_query_bin009\tAF_ref\tsame_assembly_as\n')
    for i, r in enumerate(sk[:20], 1):
        a = r['acc']; m = new.get(a, {}); twin = [b for b in byc[chk[a]] if b != a]
        f.write(f"{i}\t{a}\t{old[a]['Organism Name']}\t{old[a]['Organism Infraspecific Names Strain'] or old[a]['Organism Infraspecific Names Isolate']}\t{old[a]['Assembly Level']}\t{m.get('Assembly Release Date','')}\t{r['ANI']}\t{r['Align_fraction_query']}\t{r['Align_fraction_ref']}\t{','.join(twin)}\n")
print(open(f'{V}/top20.tsv').read())
# (vi) range 89.45-89.75
inr = [r['acc'] for r in sk if 89.45 <= float(r['ANI']) <= 89.75]
print(f'ANI 89.45-89.75: {len(inr)} records = {len({uid[a] for a in inr})} unique assemblies')
print('records 13 and 14 by rank:', [(r['acc'], r['ANI']) for r in sk[12:15]])
# full record list (i)
with open(f'{V}/records_67.tsv', 'w') as f:
    f.write('accession\torganism\tstrain_isolate\tlevel\trelease_date\tsource\tpaired_accession\tcontent_md5\tskani_ANI\n')
    ani = {r['acc']: r['ANI'] for r in sk}
    for a in sorted(old):
        m = new.get(a, {})
        f.write(f"{a}\t{old[a]['Organism Name']}\t{old[a]['Organism Infraspecific Names Strain'] or old[a]['Organism Infraspecific Names Isolate']}\t{old[a]['Assembly Level']}\t{m.get('Assembly Release Date','')}\t{old[a]['Source Database']}\t{m.get('Assembly Paired Assembly Accession','')}\t{chk[a]}\t{ani.get(a,'not reported')}\n")
for a in sorted(set(new) - set(old)):
    print('NEW:', a, new[a]['Organism Name'], new[a]['Organism Infraspecific Names Strain'] or new[a]['Organism Infraspecific Names Isolate'], new[a]['Assembly Level'], new[a]['Assembly Release Date'], 'paired:', new[a]['Assembly Paired Assembly Accession'])
