#!/usr/bin/env python3
"""Secao 17 do README (22/10/2026): rRNA/tRNA dos 25 bins e reclassificacao MIMAG. Nao altera resultados existentes.

Entradas: 16_rrna_trna/raw/<bin>/ (de scripts/rrna_trna_16.sh), DRAM annotation/{rrnas,trnas}.tsv,
          results/tables_text/TableS1_bins.tsv, FASTA dos bins (DAS Tool).
Saidas (16_rrna_trna/): registered_bin009_bin021.tsv, rrna_rerun_all.tsv, trna_rerun_all.tsv,
          compare_vs_TableS1.tsv, mimag_old_vs_new.tsv, SUMMARY.txt

Definicoes:
  cobertura do modelo = (mdl_to - mdl_from + 1) / CLEN do modelo Rfam usado pelo barrnap 1.10 (cmsearch --rfam -g, E <= 0.001);
  parcial ("p") = cobertura < 80% (mesmo limiar das notas "aligned only X percent" do barrnap dentro do DRAM; versao do barrnap
  embutida no app DRAM do KBase nao registrada, app executado em 28/08/2026; README secao 31);
  barrnap local: pacote bioconda 1.10.6, executavel se apresenta como 1.10.5;
  'trunc' do cmsearch e distancia a borda do contig registrados a parte.
  tipos de tRNA = aminoacidos distintos, excluindo Undet, SeC, Sup, fMet, iMet, Ile2 e pseudogenes (regra da Table S1 + pseudo).
  MIMAG (Bowers et al. 2017): HQ = compl >= 90, contam < 5, 5S+16S+23S presentes e >= 18 tipos de tRNA; MQ = compl >= 50 e
  contam < 10; LQ = demais. "estrita" (pedido da autora) = HQ exige 5S, 16S e 23S nao parciais.
"""
import os, re, csv, glob, collections
ROOT = os.path.expanduser('~/sbria_anammox_shotgun'); RC = ROOT + '/results/review_checks'
OUT = RC + '/16_rrna_trna/'; RAW = OUT + 'raw/'
BINS = ROOT + '/das_tool_result (1)/das_tool_output_dir_DASTool_bins/'
CMDB = os.path.expanduser('~/miniconda3/envs/rrna/db/')
EXCL = {'Undet', 'SeC', 'Sup', 'fMet', 'iMet', 'Ile2', 'Pseudo'}
GENE = {'5_8S_rRNA': '5.8S (modelo arc; fora da MIMAG)', '5S_rRNA': '5S', 'SSU_rRNA_bacteria': '16S', 'LSU_rRNA_bacteria': '23S', 'SSU_rRNA_archaea': '16S', 'LSU_rRNA_archaea': '23S'}
rep = []
def say(s=''): rep.append(s); print(s)
def w(fn, hdr, rows):
    with open(OUT + fn, 'w') as f:
        f.write('\t'.join(hdr) + '\n')
        for r in rows: f.write('\t'.join(str(x) for x in r) + '\n')

def clen(king):
    d, name = {}, None
    for l in open(f'{CMDB}{king}/{king}.rRNA.cm'):
        if l.startswith('NAME'): name = l.split()[1]
        elif l.startswith('CLEN') and name: d[name] = int(l.split()[1]); name = None
    return d
CLEN = {'bac': clen('bac'), 'arc': clen('arc')}
def contig_lens(b):
    L, k = {}, None
    for l in open(BINS + b + '.fasta'):
        if l[0] == '>': k = l[1:].split()[0]; L[k] = 0
        else: L[k] += len(l.strip())
    return L
node = lambda x: re.match(r'(NODE_\d+)', x).group(1)
bins = sorted(os.path.basename(p)[:-6] for p in glob.glob(BINS + 'bin.*.fasta'))
ts1 = {r['Bin']: r for r in csv.DictReader(open(ROOT + '/results/tables_text/TableS1_bins.tsv'), delimiter='\t')}
say('Versoes: ' + open(OUT + 'versions.txt').read().replace('\n', '; '))

# ---------- 1/2. valores registrados (DRAM) para bin.009 e bin.021
say('== 1/2. Registrado no DRAM (barrnap 0.9 e tRNAscan-SE do DRAM/KBase)')
reg = []
for b in ('bin.009', 'bin.021'):
    L = contig_lens(b)
    for r in csv.DictReader(open(ROOT + '/DRAM annotation/rrnas.tsv'), delimiter='\t'):
        if not r['fasta'].startswith(b + '.'): continue
        c = r['scaffold'].split('_assembly_')[1]; s, e = int(r['begin']), int(r['end']); cl = L[c]
        m = re.search(r'aligned only (\d+) percent', r['note'])
        reg.append((b, r['type'], node(c), cl, s, e, e - s + 1, m.group(1) if m else '>=80 (sem nota)', 'yes' if m else 'no', s - 1, cl - e))
        say(f'  {b} {r["type"]}: {node(c)} ({cl} pb) {s}-{e} ({e - s + 1} pb); nota: {r["note"] or "-"}; distancia as bordas {s - 1} / {cl - e} pb')
    tr = [r for r in csv.DictReader(open(ROOT + '/DRAM annotation/trnas.tsv'), delimiter='\t') if r['fasta'].startswith(b + '.')]
    types = {r['Type'].strip() for r in tr if r['Type'].strip() not in EXCL and 'pseudo' not in r['Note']}
    other = collections.Counter(r['Type'].strip() for r in tr if r['Type'].strip() in EXCL)
    pseudo = [f"{r['Type'].strip()} ({r['Codon']})" for r in tr if 'pseudo' in r['Note']]
    sup = [r['Codon'] for r in tr if r['Type'].strip() == 'Sup']
    say(f'  {b} tRNA: {len(tr)} genes; {len(types)} tipos de aminoacido; excluidos por tipo: {dict(other)}; pseudogenes: {pseudo or "nenhum"};'
        f' SeC: {other.get("SeC", 0)}; Sup anticodon: {sup}')
w('registered_bin009_bin021.tsv', ['bin', 'rRNA', 'contig', 'contig_len', 'begin', 'end', 'gene_len', 'barrnap09_pct_aligned', 'partial_note',
                                    'dist_start', 'dist_end'], reg)
say()

# ---------- 3. barrnap 1.10 + tRNAscan-SE 2.0.13 em todos os contigs
rr_rows, tr_rows, newsum = [], [], {}
for b in bins:
    d = RAW + b + '/'; king = 'arc' if b == 'bin.003' else 'bac'; L = contig_lens(b)
    tbl = [l.split() for f in glob.glob(d + '*.rRNA.tbl') for l in open(f) if not l.startswith('#')]
    calls = [l.split('\t') for l in open(d + 'barrnap.gff') if l.strip() and not l.startswith('#')]
    best = {}
    for g in calls:
        nm = re.search(r'Alias=([^;]+)', g[8]).group(1); s, e = int(g[3]), int(g[4])
        hit = [t for t in tbl if t[0] == g[0] and t[2] == nm and min(int(t[7]), int(t[8])) == s and max(int(t[7]), int(t[8])) == e]
        t = hit[0] if hit else None
        cov = 100 * (int(t[6]) - int(t[5]) + 1) / CLEN[king][nm] if t else float('nan')
        gene = GENE[nm]; cl = L[g[0]]
        row = (b, gene, node(g[0]), cl, s, e, e - s + 1, CLEN[king][nm], f'{cov:.1f}', t[10] if t else '?', 'p' if cov < 80 else 'complete',
               s - 1, cl - e, t[15] if t else '')
        rr_rows.append(row)
        if gene not in best or cov > float(best[gene][8]): best[gene] = row
    tr = []
    if os.path.exists(d + 'trnascan.txt'):
        for l in open(d + 'trnascan.txt').readlines()[3:]:
            p = [x.strip() for x in l.rstrip('\n').split('\t')]
            if len(p) >= 9: tr.append(p)
    for p in tr: tr_rows.append([b] + p)
    types = {p[4] for p in tr if p[4] not in EXCL and 'pseudo' not in (p[9] if len(p) > 9 else '')}
    newsum[b] = dict(r={k: best[k][10] for k in best}, cov={k: best[k][8] for k in best}, n=len(tr), types=len(types),
                     sec=sum(p[4] == 'SeC' for p in tr), sup=[p[5] for p in tr if p[4] == 'Sup'],
                     pseudo=[f'{p[4]} ({p[5]})' for p in tr if len(p) > 9 and 'pseudo' in p[9]])
w('rrna_rerun_all.tsv', ['bin', 'rRNA', 'contig', 'contig_len', 'begin', 'end', 'gene_len', 'model_len', 'model_cov_pct', 'cmsearch_trunc',
                         'status', 'dist_start', 'dist_end', 'evalue'], rr_rows)
w('trna_rerun_all.tsv', ['bin', 'contig', 'tRNA_n', 'begin', 'end', 'type', 'anticodon', 'intron_begin', 'intron_end', 'score', 'note'], tr_rows)

def fmt(b):
    s = newsum[b]; out = []
    for k in ('5S', '16S', '23S'): out.append(k + (' p' if s['r'].get(k) == 'p' else '') if k in s['r'] else '–')
    return '/'.join(out)
cmp_rows, changed = [], []
for b in bins:
    old_r, old_t = ts1[b]['rRNA genes (5S/16S/23S)'], int(ts1[b]['tRNA types']); nr, nt = fmt(b), newsum[b]['types']
    ch = (old_r.replace(' ', '') != nr.replace(' ', '')) or old_t != nt
    cmp_rows.append((b, old_r, nr, old_t, nt, newsum[b]['n'], 'yes' if ch else 'no'))
    if ch: changed.append(f'{b}: rRNA {old_r} -> {nr}; tRNA types {old_t} -> {nt}')
w('compare_vs_TableS1.tsv', ['bin', 'TableS1_rRNA', 'rerun_rRNA', 'TableS1_tRNA_types', 'rerun_tRNA_types', 'rerun_tRNA_genes', 'changed'], cmp_rows)
say('== 3. barrnap 1.10.5 + tRNAscan-SE 2.0.13 em todos os contigs; bins com mudanca em relacao a Table S1:')
for c in changed: say('  ' + c)
if not changed: say('  nenhum')
say('Para bin.009 e bin.021 (rerun):')
for r in rr_rows:
    if r[0] in ('bin.009', 'bin.021'): say(f'  {r[0]} {r[1]}: {r[2]} ({r[3]} pb) {r[4]}-{r[5]} ({r[6]} pb); modelo {r[7]} -> {r[8]}%; trunc={r[9]}; {r[10]}; bordas {r[11]}/{r[12]} pb')
for b in ('bin.009', 'bin.021'):
    s = newsum[b]; say(f'  {b} tRNA: {s["n"]} genes; {s["types"]} tipos; SeC {s["sec"]}; Sup {s["sup"]}; pseudo {s["pseudo"] or "nenhum"}')
say()

# ---------- 4. MIMAG antiga x nova
def tier(b, strict):
    c, k = float(ts1[b]['Completeness (%)']), float(ts1[b]['Contamination (%)']); s = newsum[b]
    rr = all(g in s['r'] for g in ('5S', '16S', '23S')) and (not strict or all(s['r'][g] == 'complete' for g in ('5S', '16S', '23S')))
    if c >= 90 and k < 5 and rr and s['types'] >= 18: return 'High'
    return 'Medium' if c >= 50 and k < 10 else 'Low'
rows = []
say('== 4. MIMAG: Table S1 x rerun (regra publicada) x rerun (estrita: rRNA nao parciais)')
for b in bins:
    old, pub, st = ts1[b]['Tier (MIMAG)'], tier(b, False), tier(b, True)
    rows.append((b, ts1[b]['Completeness (%)'], ts1[b]['Contamination (%)'], ts1[b]['rRNA genes (5S/16S/23S)'], fmt(b), ts1[b]['tRNA types'],
                 newsum[b]['types'], old, pub, st))
    hq = float(ts1[b]['Completeness (%)']) >= 90 and float(ts1[b]['Contamination (%)']) < 5
    if hq or old != pub or old != st:
        say(f'  {b} ({ts1[b]["Completeness (%)"]}% / {ts1[b]["Contamination (%)"]}%): rRNA {fmt(b)}, {newsum[b]["types"]} tipos de tRNA -> '
            f'antiga {old}; publicada {pub}; estrita {st}' + ('  <-- muda' if old != pub or old != st else ''))
w('mimag_old_vs_new.tsv', ['bin', 'completeness', 'contamination', 'TableS1_rRNA', 'rerun_rRNA', 'TableS1_tRNA_types', 'rerun_tRNA_types',
                           'tier_TableS1', 'tier_rerun_published', 'tier_rerun_strict'], rows)
for col, nm in ((7, 'Table S1'), (8, 'rerun publicada'), (9, 'rerun estrita')):
    say(f'  contagem {nm}: {dict(collections.Counter(r[col] for r in rows))}')
open(OUT + 'SUMMARY.txt', 'w').write('\n'.join(rep) + '\n')
