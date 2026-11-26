# Valores da Figure S2 — versão 19/10/2026 (anterior: make_FigS2_values_v_antiga2.py).
# Fonte única: Table S2 (results/review_checks/supp_tables/TableS2_targeted_hits.tsv): identidade, cobertura da
# referência e status pelo critério dos Methods (e-value <= 1e-10; identidade >= 40%; cobertura >= 70%; parálogos
# excluídos). Comprimento da referência: aba References (TableS2_references.tsv).
# Status DRAM/KOfam do gene do hit (gene DRAM sobreposto ao locus, coluna "DRAM gene at locus" da S2):
#   Yes     = gene com o KO esperado
#   Partial = gene só com domínio Pfam PF00384 (molibdopterina), sem KO
#   No      = sem anotação correspondente, ou contig fora de bin consolidado
# Barras novas (19/10/2026): amoA/amoB de NODE_13511 (atribuídos ao bin.007, ‡), amoC de NODE_2578 (bin.007),
# nosZ do bin.017; amoC de NODE_4969 passa a usar a AmoC completa de N. europaea (WP_011111554.1, 271 aa).
import os, csv, collections

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
os.chdir(ROOT)
SUP = 'results/review_checks/supp_tables/'
S2 = list(csv.DictReader(open(SUP + 'TableS2_targeted_hits.tsv'), delimiter='\t'))
REFLEN = {}
for r in csv.DictReader(open(SUP + 'TableS2_references.tsv'), delimiter='\t'):
    if r['Length (aa)']: REFLEN.setdefault(r['Accession'], int(float(r['Length (aa)'])))
PFAM = {}
for r in csv.reader(open('DRAM annotation/annotations.tsv'), delimiter='\t'):
    if r[0] and r[1] != 'fasta':
        n = r[2].split('_length')[0]; PFAM[n + '_' + r[0].rsplit('_', 1)[1]] = (r[8], r[17])

# (gene S2, rótulo do eixo, organismo, seleção: lista de contigs ou bin, KO esperado, nota extra)
spec = [
    ('amoA', 'amoA', 'Nitrosomonas (NODE_4969)', ['NODE_4969'], 'K10944', ''),
    ('amoC', 'amoC', 'Nitrosomonas (NODE_4969)', ['NODE_4969'], 'K10946', ''),
    ('amoA', 'amoA (comammox)', 'Nitrospira_F (bin.007‡)', ['NODE_13511'], 'K10944', 'indirect'),
    ('amoB', 'amoB (comammox)', 'Nitrospira_F (bin.007‡)', ['NODE_13511'], 'K10945', 'indirect'),
    ('amoC', 'amoC (comammox)', 'Nitrospira_F (bin.007)', ['NODE_2578'], 'K10946', ''),
    ('nxrA', 'nxrA', 'Nitrospira_F (bin.007)', ['NODE_4141', 'NODE_348'], 'K00370', ''),
    ('nxrA', 'nxrA (fragment, 135 aa)', 'Nitrospira_A (bin.004)', ['NODE_719'], 'K00370', ''),
    ('nxrA', 'nxrA (unbinned)', 'Nitrospira (NODE_2432)', ['NODE_2432'], 'K00370', ''),
    ('hzsA', 'hzsA', 'Anammox contig (NODE_5668)', ['NODE_5668'], None, ''),
    ('hzsB', 'hzsB', 'Anammox contig (NODE_5312)', ['NODE_5312'], None, ''),
    ('hzsC', 'hzsC', 'Anammox contig (NODE_7170)', ['NODE_7170'], None, ''),
    ('hdh', 'hdh', 'Anammox contig (NODE_23806)', ['NODE_23806'], None, ''),
    ('nirS', 'nirS', 'Ca. Brocadia (bin.009)', ['NODE_696'], 'K15864', ''),
    ('nxrA', 'nxrA (anammox)', 'Ca. Brocadia (bin.009)', ['NODE_4298', 'NODE_1742'], 'K00370', ''),
    ('narG', 'narG', 'VBCG01 (bin.008)', 'bin.008', 'K00370', 'split gene: 2 adjacent ORFs on NODE_2924 (DRAM)'),
    ('narG', 'narG', 'Ca. Desulfobacillus (bin.011)', 'bin.011', 'K00370', ''),
    ('napA', 'napA', 'VBCG01 (bin.008)', 'bin.008', 'K02567', ''),
    ('napA', 'napA', 'Sphaerotilus (bin.016)', 'bin.016', 'K02567', ''),
    ('napA', 'napA', 'Ca. Desulfobacillus (bin.011)', 'bin.011', 'K02567', ''),
    ('napA', 'napA', 'Thiothrix (bin.017)', 'bin.017', 'K02567', ''),
    ('nosZ', 'nosZ (clade II)', 'Ca. Desulfobacillus (bin.011)', 'bin.011', 'K00376', ''),
    ('nosZ', 'nosZ (clade II)', 'Thiothrix (bin.017)', 'bin.017', 'K00376', ''),
]
RANK = {'supported': 0, 'partial': 1}
f = lambda x: float(x) if x not in ('', None) else float('nan')
rows = []
for gene, axis, org, sel, ko, extra in spec:
    cand = [r for r in S2 if r['Gene'] == gene and r['Evidence'] in RANK and
            ((r['Contig'] in sel) if isinstance(sel, list) else r['Bin'].rstrip('†‡') == sel)]
    b = min(cand, key=lambda r: (RANK[r['Evidence']], -f(r['Bit score'])))
    ref = b['Reference (accession)']; L = REFLEN[ref]; cov = f(b['Reference coverage (%)'])
    others = sorted({(r['Contig'], f(r['Identity (%)'])) for r in cand if r is not b and r['Evidence'] == 'supported'
                     and r['Contig'] != b['Contig']}, key=lambda x: -x[1])
    ko_, pf = PFAM.get(b['DRAM gene at locus'], ('', ''))
    dram = 'Yes' if ko and ko_ == ko else ('Partial' if (not ko_ and 'PF00384' in pf) else 'No')
    note = '; '.join(x for x in [extra, b['Note']] if x)
    rows.append([gene, axis, org, f"{f(b['Identity (%)']):.2f}", str(round(cov * L / 100)), str(L), f"{cov:.1f}",
                 'TRUE' if b['Evidence'] == 'partial' else 'FALSE', b['Evidence'], b['E-value'], b['Contig'], ref,
                 'TableS2_targeted_hits.tsv', dram, '; '.join(f'{n} {p:.1f}%' for n, p in others), note])

hdr = ['gene', 'axis_label', 'organism', 'identity_pct', 'aln_len', 'ref_len', 'ref_coverage_pct', 'partial', 'blast_status',
       'evalue', 'contig', 'reference', 'source_file:line', 'dram', 'other_full_length_copies', 'note']
with open('results/tables_text/FigS2_values.tsv', 'w') as o:
    o.write('\t'.join(hdr) + '\n')
    for r in rows: o.write('\t'.join(r) + '\n')
for r in rows: print('\t'.join(r[:3] + [r[3], r[6], r[8], r[11], r[13], r[14]]))
