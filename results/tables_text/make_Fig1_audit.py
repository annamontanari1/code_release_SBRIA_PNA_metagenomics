# Auditoria de evidência da Figure 1 (organismo x gene) — versão 19/10/2026.
# Fonte única: Table S2 (results/review_checks/supp_tables/TableS2_targeted_hits.tsv, gerada por make_TableS2.py),
# para que Figura 1, Figura S2 e Table S2 mostrem os mesmos valores. Versão anterior (recalculava a partir de
# results/targeted_search/ com AmoB/AmoC truncadas e só NosZ de clado I): make_Fig1_audit_v_antiga2.py.
# Critério (Methods): e-value <= 1e-10; identidade >= 40%; cobertura da referência >= 70% -> supported, < 70% -> partial;
# parálogos (hdh/hao, napA/fdh, narG/nxrA) excluídos por KO discordante ou, para hdh, por melhor hit recíproco HAO em
# K. stuttgartiensis (aba "Excluded paralog hits" da Table S2).
# Níveis por célula: supported / partial / kofam (gene-alvo anotado pelo DRAM sem acerto dirigido aceito) /
# not detected / not tested (só o narG/nxrA family do bin.009, mostrado na coluna nxrA).
# Unidades não binadas (NODE_4969, contigs anammox, NODE_2432): genes de desnitrificação (nirS, narG/nxrA family, napA,
# norB, norC, nosZ) ficam "not tested", como nas versões anteriores da figura (a busca original desses genes foi feita
# contra os genes dos bins); pedido da autora, 19/10/2026 (versão anterior: make_Fig1_audit_v_antiga3.py).
# Linha nova: Nitrospira (NODE_2432, unbinned), nxrA.
# bin.007: amoA/amoB de NODE_13511 (não binado) atribuídos ao bin.007 (profundidade, GC, TNF, co-binning MaxBin2 e
# 46 pares de reads; results/review_checks/comammox/) — marcados como atribuição indireta (‡ na figura).
import os, csv, collections

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
os.chdir(ROOT)
S2 = list(csv.DictReader(open('results/review_checks/supp_tables/TableS2_targeted_hits.tsv'), delimiter='\t'))

# KOfam (DRAM) por bin
KOS = collections.defaultdict(set)
for r in csv.reader(open('DRAM annotation/annotations.tsv'), delimiter='\t'):
    if r[0] and r[1] != 'fasta' and r[8]: KOS[r[1].split('.fasta')[0]].add(r[8])

ANAMMOX = {'NODE_5668', 'NODE_5312', 'NODE_7170', 'NODE_17695', 'NODE_29475', 'NODE_47442', 'NODE_20914', 'NODE_24842', 'NODE_23806'}
units = {
    'Nitrosomonas (NODE_4969)': dict(contigs={'NODE_4969'}),
    'Nitrospira_A (bin.004)': dict(bin='bin.004'),
    'Nitrospira_F (bin.007)': dict(bin='bin.007', extra={'NODE_13511'}),
    'Nitrospira (NODE_2432, unbinned)': dict(contigs={'NODE_2432'}),
    'Ca. Brocadia (bin.009)': dict(bin='bin.009'),
    'Anammox contigs (unassigned)': dict(contigs=ANAMMOX),
    'bin.011': dict(bin='bin.011'), 'bin.016': dict(bin='bin.016'),
    'bin.008': dict(bin='bin.008'), 'bin.017': dict(bin='bin.017'),
}
NXR_UNITS = {'Nitrospira_A (bin.004)', 'Nitrospira_F (bin.007)', 'Nitrospira (NODE_2432, unbinned)', 'Ca. Brocadia (bin.009)'}
DENITRIF = {'nirS', 'narG/nxrA family', 'napA', 'norB', 'norC', 'nosZ'}
GENES = [('amoA', 'K10944'), ('amoB', 'K10945'), ('amoC', 'K10946'), ('nxrA', 'K00370'), ('hzsA', None), ('hzsB', None),
         ('hzsC', None), ('hdh', None), ('nirS', 'K15864'), ('narG/nxrA family', 'K00370'), ('napA', 'K02567'),
         ('norB', 'K04561'), ('norC', 'K02305'), ('nosZ', 'K00376')]
S2GENE = {'nxrA': {'nxrA'}, 'narG/nxrA family': {'narG', 'nxrA', 'narG/nxrA (unassigned)'}}
RANK = {'supported': 0, 'partial': 1, 'KOfam only': 2}
EXTRA_NOTE = {('Ca. Brocadia (bin.009)', 'nxrA'): 'anammox-type NXR; reference Ca. K. stuttgartiensis NxrA (kustd1700); Nitrospira NxrA 58.3%'}

def fnum(x):
    try: return float(x)
    except (TypeError, ValueError): return float('nan')

out = []
for u, ud in units.items():
    for gene, ko in GENES:
        if (gene == 'nxrA' and u not in NXR_UNITS) or (gene == 'narG/nxrA family' and u in NXR_UNITS): continue
        names = S2GENE.get(gene, {gene}); binned = 'bin' in ud
        def in_unit(r):
            if 'bin' in ud:
                return r['Bin'].rstrip('†‡') == ud['bin'] or r['Contig'] in ud.get('extra', set())
            return r['Contig'] in ud['contigs']
        cand = [r for r in S2 if r['Gene'] in names and r['Evidence'] in RANK and in_unit(r)]
        best = min(cand, key=lambda r: (RANK[r['Evidence']], -fnum(r['Bit score']))) if cand else None
        kofam = binned and ko in KOS[ud['bin']]
        if not binned and gene in DENITRIF: ev = 'not tested'
        elif best and best['Evidence'] in ('supported', 'partial'): ev = best['Evidence']
        elif kofam or (best and best['Evidence'] == 'KOfam only'): ev = 'kofam'
        else: ev = 'not detected'
        b = best if ev in ('supported', 'partial') else None
        note = b['Note'] if b else (best['Note'] if best else '')
        if b and b['Bin'].endswith(('†', '‡')): note = 'indirect attribution to bin.007; ' + note
        note = '; '.join(x for x in [note, EXTRA_NOTE.get((u, gene), '')] if x)
        out.append(dict(unit=u, gene=gene, evidence=ev,
                        blast_status=b['Evidence'] if b else ('KOfam only' if ev == 'kofam' else ('not searched' if ev == 'not tested' else 'no hit')),
                        best_hit_contig=b['Contig'] if b else '', identity_pct=b['Identity (%)'] if b else '',
                        evalue=b['E-value'] if b else '', ref_coverage_pct=b['Reference coverage (%)'] if b else '',
                        reference=b['Reference (accession)'] if b else '', source='TableS2_targeted_hits.tsv' if b else '',
                        KOfam=(ko or 'n.a.') + (':yes' if kofam else (':no' if binned else ':n.a.')),
                        DRAM_annotation_of_hit_gene=(b['DRAM KO'] if b else ''), note=note))

out.append({**{k: '' for k in out[0]}, 'unit': 'Ca. Brocadia (bin.009)', 'gene': 'narG/nxrA family', 'evidence': 'not tested',
            'blast_status': 'not applicable', 'KOfam': 'K00370:n.a.',
            'note': 'not applicable: NxrA/NarG-family genes of bin.009 shown under nxrA'})
out.append({**{k: '' for k in out[0]}, 'unit': 'Nitrospira (NODE_2432, unbinned)', 'gene': 'narG/nxrA family',
            'evidence': 'not tested', 'blast_status': 'not applicable', 'KOfam': 'K00370:n.a.',
            'note': 'not applicable: unbinned contig; NxrA of NODE_2432 shown under nxrA'})
cols = list(out[0].keys())
with open('results/tables_text/Fig1_evidence_audit.tsv', 'w') as o:
    o.write('\t'.join(cols) + '\n')
    for r in out: o.write('\t'.join(r[c] for c in cols) + '\n')

# ---- diferenças em relação à versão anterior ----------------------------------
old = {(r['unit'], r['gene']): r for r in
       csv.DictReader(open('results/tables_text/Fig1_evidence_audit_v_antiga3.tsv'), delimiter='\t')}
print('CHANGED CELLS (old -> new):')
for r in out:
    o_ = old.get((r['unit'], r['gene']))
    ov = o_['evidence'] if o_ else '(absent)'
    so = o_['identity_pct'] if o_ and ov in ('supported', 'partial') else ''
    sn = r['identity_pct'] if r['evidence'] in ('supported', 'partial') else ''
    if ov != r['evidence'] or (so and sn and abs(float(so) - float(sn)) >= 0.05) or bool(so) != bool(sn):
        print(f"{r['unit']}\t{r['gene']}\t{ov} {so}\t->\t{r['evidence']} {sn}\t{r['best_hit_contig']}\t{r['note'][:80]}")
