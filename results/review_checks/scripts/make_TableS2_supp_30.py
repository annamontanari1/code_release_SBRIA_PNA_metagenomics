#!/usr/bin/env python3
"""Seção 30 do README (24/10/2026): Table S2 do material suplementar -> results/supplementary_tables/ (Table_S2.xlsx + um CSV por aba).
Só lê arquivos existentes; não recalcula identidades nem coberturas (a aba S2D usa a saída de 29_panel_masking0/, seção 29).
Fontes: supp_tables/TableS2_targeted_hits.tsv e TableS2_references.tsv (busca dirigida, 21/10/2026), 14_ncycdb/TableS2_vs_NCycDB.tsv,
19_hzsA_v2/hzsA_reference_candidates.*, 12_anammox_contigs/anammox_reference_panel.tsv, 11_sulfur/dsrA_tree_input.faa,
06_tableS2/hzs_orfs_vs_anammox_panel_long.tsv, DRAM annotation/annotations.tsv.
Uso (de ~/sbria_anammox_shotgun): python results/review_checks/scripts/make_TableS2_supp_30.py
"""
import os, re, glob
import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Font, Alignment

RC = 'results/review_checks'
OUT = 'results/supplementary_tables'
os.makedirs(OUT, exist_ok=True)
r1 = lambda x: '' if pd.isna(x) or x == '' else f'{float(x):.1f}'
FIX = {'KKO20885.1': 'hzsB', 'KKO20886.1': 'hzsC'}

# ------------------------------------------------------------------ cabeçalhos FASTA (organismo do registro NCBI)
HDR = {}
for f in [x for x in glob.glob(f'{RC}/refs/*.faa') if 'okubo_nosZ' not in x] + [f'{RC}/12_anammox_contigs/anammox_refs.faa', f'{RC}/19_hzsA_v2/hzsA_reference_candidates.faa',
                                           f'{RC}/refs/dsrA_refs_panel.faa']:
    for l in open(f):
        if not l.startswith('>'):
            continue
        h = l[1:].strip()
        m = re.match(r'(?:sp|tr)\|([^|]+)\|', h)
        acc = m.group(1) if m else h.split()[0].split('|')[-1] if h.split() else ''
        org = re.findall(r'\[([^\[\]]+)\]', h.split(' >')[0])
        if acc and org and not org[-1].startswith(('locus_tag', 'protein')):
            HDR.setdefault(acc, org[-1])
for _, x in pd.read_csv(f'{RC}/12_anammox_contigs/anammox_reference_panel.tsv', sep='\t').iterrows():
    o = re.findall(r'\[([^\[\]]+)\]', x.header)
    if o:
        HDR.setdefault(x.accession, o[-1])
MULTI = set()
for f in glob.glob(f'{RC}/refs/*.faa') + [f'{RC}/19_hzsA_v2/hzsA_reference_candidates.faa']:
    MULTI |= {l[1:].split()[0] for l in open(f) if l.startswith('>') and 'MULTISPECIES' in l}
MULTI |= set(pd.read_csv(f'{RC}/12_anammox_contigs/anammox_reference_panel.tsv', sep='\t').query('header.str.contains("MULTISPECIES")', engine='python').accession)
MULTI |= set(pd.read_csv(f'{RC}/supp_tables/TableS2_references.tsv', sep='\t').query('Description.str.contains("MULTISPECIES", na=False) or Organism.str.contains("MULTISPECIES", na=False)', engine='python').Accession)
HDR['BBO22106.1'] = 'Candidatus Desulfobacillus denitrificans'   # cabeçalho de okubo_nosZ.faa traz anotação livre
# 24/10/2026 (README seção 35): rifle_RCA_RCB_amo.faa traz anotação livre no cabeçalho; organismo do registro NCBI = "Nitrospira sp.", /isolate RCA ou RCB
for _a in ['TKS61416.1', 'TKS61417.1', 'TKS61418.1']: HDR[_a] = 'Nitrospira sp. RCA'
for _a in ['TKS59186.1', 'TKS59187.1', 'TKS59188.1']: HDR[_a] = 'Nitrospira sp. RCB'
SEQLEN = {}
for f in [f'{RC}/12_anammox_contigs/anammox_refs.faa', f'{RC}/19_hzsA_v2/hzsA_reference_candidates.faa']:
    k = None
    for l in open(f):
        if l.startswith('>'):
            h = l[1:].strip(); m = re.match(r'(?:sp|tr)\|([^|]+)\|', h); k = m.group(1) if m else h.split()[0]; SEQLEN[k] = 0
        else:
            SEQLEN[k] += len(l.strip())

# ------------------------------------------------------------------ S2A: References
ref = pd.read_csv(f'{RC}/supp_tables/TableS2_references.tsv', sep='\t')
blastp_set = set(l[1:].split()[0] for l in open(f'{RC}/refs/refs_all_for_blastp.faa') if l.startswith('>') and l[1:].strip())
blastp_set |= {l.strip()[1:].split('|')[1] for l in open(f'{RC}/refs/sox_Paracoccus_refs.faa') if l.startswith('>')}
blastp_set |= {'BBO20425.1', 'BBO20424.1', '8SR2_C', 'WP_256036300.1', 'WP_099325812.1'}
dsr = {}
for l in open(f'{RC}/11_sulfur/dsrA_tree_input.faa'):
    if l.startswith('>') and '|' in l:
        lab, acc = l[1:].strip().split('|')
        dsr[acc] = ('DsrA typing: test sequence, not used in the assignment' if lab.startswith('TEST') else
                    'DsrA oxidative' if lab.startswith('OX_') else 'DsrA reductive (archaeal, root)' if lab.startswith('ARCH') else 'DsrA reductive')
cand = pd.read_csv(f'{RC}/19_hzsA_v2/hzsA_reference_candidates.tsv', sep='\t')
tree12 = set(cand[cand.proposed == 'yes'].accession)
panel = pd.read_csv(f'{RC}/12_anammox_contigs/anammox_reference_panel.tsv', sep='\t')

GENE_FIX = {'YEH58028.1': 'nxrA', 'Q1PZD8': 'nxrA (anammox-type; kustd1700)', 'WP_230999932.1': 'narG', 'WP_230999923.1': 'nosZ (clade I)',
            'WP_234989093.1': 'amoC', 'KKO20885.1': 'hzsB', 'KKO20886.1': 'hzsC'}
CITED = {'WP_011111554.1': 'Nitrosomonas europaea', 'WP_099325812.1': 'Candidatus Kuenenia stuttgartiensis',
         'Q1PZD8': 'Candidatus Kuenenia stuttgartiensis', 'YEH58028.1': 'Nitrospira', 'YEH58030.1': 'Nitrospira',
         'WP_301017512.1': 'Thiothrix', 'WP_505736924.1': 'Methanothrix', 'NZ_CP049055.1': 'Candidatus Kuenenia stuttgartiensis CSTR1'}
NOTE = {**{x: 'comammox clade A (Poghosyan et al., 2019); NCBI organism "Nitrospira sp.", isolate RCA' for x in ['TKS61416.1', 'TKS61417.1', 'TKS61418.1']},
        **{x: 'comammox clade B (Poghosyan et al., 2019); NCBI organism "Nitrospira sp.", isolate RCB' for x in ['TKS59186.1', 'TKS59187.1', 'TKS59188.1']},
        'TRZ97369.1': 'organism named at family rank (NCBI "Nitrospiraceae bacterium", isolate PowLak16_MAG19)',
        'Q1PZD8': 'Kuenenia NxrA cited in the text as kustd1700 (locus tag of this UniProtKB entry); one entry, not two references',
        'KKO20885.1': 'GenBank label "hydrazine synthase subunit C" is swapped; by sequence this is HZS-beta (hzsB; PDB 5C2V chain beta, Q1Q0T4)',
        'KKO20886.1': 'GenBank label "hydrazine synthase subunit B" is swapped; by sequence this is HZS-gamma (hzsC; PDB 5C2V chain gamma, Q1Q0T3)',
        'WP_164994616.1': 'HZS-beta (Swiss-Prot/PDB 5C2V nomenclature; TrEMBL and NCycDB swap beta and gamma)',
        'WP_478210066.1': 'HZS-gamma (Swiss-Prot/PDB 5C2V nomenclature; TrEMBL and NCycDB swap beta and gamma)',
        'Q1Q0T2.1': 'query used to select the HzsA tree references (identity >=40%, coverage >=70%); identical (100%) to WP_099325812.1',
        'WP_234989093.1': 'partial AmoC record (115 aa) of the original search; replaced by WP_011111554.1 (271 aa)',
        'NZ_CP049055.1': 'genome record; its predicted proteome was used for the reciprocal best-hit test separating hdh from the HAO family',
        'BBO20425.1': 'dsrA of Ca. D. denitrificans: targeted-search reference; included in the DsrA typing tree only as a test sequence',
        'OQD44929.1': 'HAO-like protein (26-30% identity to the hdh of the contigs); listed as hdh in the panel by its GenBank label',
        'KKO18440.1': 'HAO-like protein (26-30% identity to the hdh of the contigs); listed as hdh in the panel by its GenBank label'}

OTHER_USE = {'Q82W83.1': 'other: reference set of the original amo search', 'WP_011111571.1': 'other: reference set of the original amo search',
             'P24474.1': 'other: nirS phylogeny reference', 'WP_203415591.1': 'other: nirS phylogeny reference'}

def norm_org(o):
    o = re.sub(r'\s*\((?:RefSeq MULTISPECIES|kustd1700)[^)]*\)', '', str(o)).strip()
    return o

rows = {}
def add(acc, gene, use, org=None, length=None, src=None):
    r = rows.setdefault(acc, dict(gene=GENE_FIX.get(acc, gene), uses=[], org=None, length=None, src=None))
    if use not in r['uses']:
        r['uses'].append(use)
    r['org'] = r['org'] or org; r['length'] = r['length'] or length; r['src'] = r['src'] or src

for _, x in ref.iterrows():
    a, gu = x.Accession, str(x['Gene / use'])
    gene = re.match(r'[A-Za-z/]+', gu).group(0)
    if 'clade II' in gu or a in {'BBO22106.1', 'WP_301017512.1', 'AAZ46324.1', 'AAZ46320.1', 'WP_129545366.1', 'WP_011421454.1', 'GAB1413680.1'}:
        gene = 'nosZ (clade II)'
    if gu.startswith('amo gamma'): gene = 'amo gamma subunit (cytochrome c1)'
    if '(HAO-like)' in gu: gene = 'hdh (HAO-like)'
    if 'AmoA tree' in gu:
        use = 'AmoA tree (Figure S3)'
    elif a in blastp_set:
        use = 'targeted search'
    elif a in dsr:
        use = dsr[a]
    else:
        use = OTHER_USE[a]
    org = HDR.get(a) or norm_org(x.Organism)
    add(a, gene, use, org, int(x['Length (aa)']) if pd.notna(x['Length (aa)']) else None, x.Source)
    if a in dsr and dsr[a] != use:
        add(a, gene, dsr[a])
for _, x in panel.iterrows():
    add(x.accession, FIX.get(x.accession, x.gene), 'anammox panel', HDR.get(x.accession), SEQLEN.get(x.accession),
        'RefSeq' if x.accession.startswith('WP_') else 'Swiss-Prot' if x.accession == 'Q1PW30.2' else 'GenBank (INSDC)')
for _, x in cand.iterrows():
    if x.accession in tree12 or x.accession == 'Q1Q0T2.1':
        add(x.accession, 'hzsA', 'HzsA tree' if x.accession in tree12 else 'HzsA tree: reference-selection query',
            HDR.get(x.accession) or ('Candidatus Kuenenia stuttgartiensis' if x.accession == 'Q1Q0T2.1' else None), int(x.length_aa),
            {'RefSeq': 'RefSeq', 'GenBank': 'GenBank (INSDC)', 'Swiss-Prot': 'Swiss-Prot'}[x.source])
add('NZ_CP049055.1', 'hdh/HAO family (proteome)', 'other: proteome for the hdh/HAO reciprocal best-hit test',
    'Candidatus Kuenenia stuttgartiensis', None, 'RefSeq (nucleotide, complete genome)')
for acc in ['WP_024845874.1', 'WP_024845875.1', 'WP_024845876.1', 'WP_024845877.1', 'WP_024845878.1', 'WP_024845879.1', 'WP_024845880.1']:
    rows[acc]['org'] = 'Paracoccus'   # registros RefSeq MULTISPECIES [Paracoccus]
rows['Q1PZD8']['org'] = 'Candidatus Kuenenia stuttgartiensis'
rows['8SR2_C']['src'] = 'PDB'

USE_ORDER = ['targeted search', 'anammox panel', 'HzsA tree', 'HzsA tree: reference-selection query', 'DsrA oxidative', 'DsrA reductive',
             'DsrA reductive (archaeal, root)', 'DsrA typing: test sequence, not used in the assignment', 'AmoA tree (Figure S3)']
S2A = []
for a, r in rows.items():
    rec = r['org'] or ''
    cited = CITED.get(a, rec.replace('stuttgartensis', 'stuttgartiensis'))
    note = NOTE.get(a, '')
    if a in MULTI:
        note = '; '.join(x for x in ['RefSeq MULTISPECIES record (organism given at the rank of the record)', note] if x)
    if 'stuttgartensis' in rec:
        note = '; '.join(x for x in [note, 'NCBI spelling "stuttgartensis"'] if x)
    uses = sorted(r['uses'], key=lambda u: USE_ORDER.index(u) if u in USE_ORDER else 99)
    S2A.append({'Gene/protein': r['gene'], 'NCBI record organism': rec, 'Reference organism cited': cited, 'Accession': a,
                'Use': '; '.join(uses), 'Source database': r['src'] or '', 'Length (aa)': '' if r['length'] is None else int(r['length']),
                'Note': note})
S2A = pd.DataFrame(S2A).sort_values(['Gene/protein', 'Reference organism cited', 'Accession'], key=lambda s: s.str.lower()).reset_index(drop=True)
assert S2A.Accession.is_unique

# ------------------------------------------------------------------ S2B: Targeted search and NCycDB
t = pd.read_csv(f'{RC}/supp_tables/TableS2_targeted_hits.tsv', sep='\t', dtype=str).fillna('')
n = pd.read_csv(f'{RC}/14_ncycdb/TableS2_vs_NCycDB.tsv', sep='\t', dtype=str).fillna('')
assert len(t) == len(n) == 392 and (t[['Gene', 'Contig', 'Evidence']].values == n[['Gene', 'Contig', 'Evidence']].values).all()
KO = {'amoA': 'K10944', 'amoB': 'K10945', 'amoC': 'K10946', 'narG': 'K00370', 'nxrA': 'K00370', 'narG/nxrA (unassigned)': 'K00370',
      'nxrA/narG': 'K00370', 'nxrB': 'K00371', 'nxrB-like': 'K00371', 'hdh': 'K20935', 'napA': 'K02567', 'nirS': 'K15864', 'nirK': 'K00368',
      'nrfA': 'K03385', 'nirB': 'K00362', 'nirD': 'K00363', 'norB': 'K04561', 'norC': 'K02305', 'nosZ': 'K00376', 'mcrA': 'K00399',
      'soxA': 'K17222', 'soxB': 'K17224', 'soxC': 'K17225', 'soxD': 'K22622', 'soxX': 'K17223', 'soxY': 'K17226', 'soxZ': 'K17227',
      'dsrA': 'K11180', 'dsrB': 'K11181'}
PFAM = {}
for l in open('DRAM annotation/annotations.tsv'):
    r = l.rstrip('\n').split('\t')
    if r[0] and r[1] != 'fasta':
        PFAM[r[2].split('_length')[0] + '_' + r[0].rsplit('_', 1)[1]] = (r[8], r[17])
def dram(row):
    if row['Evidence'] == 'not detected' or not row['DRAM gene at locus']:
        return 'no'
    ko, pf = PFAM.get(row['DRAM gene at locus'], (row['DRAM KO'], ''))
    if ko and ko == KO.get(row['Gene']):
        return 'yes'
    return 'partial' if (not ko and 'PF00384' in pf) else 'no'
OUT_SCOPE = re.compile(r'^(sox|dsr|mcr)')
def status(i):
    if t.Evidence[i] != 'supported':
        return 'n.a.'
    if OUT_SCOPE.match(t.Gene[i]):
        return 'out of scope'
    res, ex = n.result[i], n.discordance_explanation[i]
    return {'concordant': 'concordant', 'same_gene_below_threshold': 'same gene below threshold',
            'homolog_not_in_map': 'same-family homolog outside the map', 'other_gene_below_threshold': 'other gene below threshold'}.get(res) or \
        {'pmo': 'pmo (amo-assigned)', 'nar/nxr family': 'nar/nxr family', 'hzsB/hzsC nomenclature': 'hzsB/hzsC nomenclature swap'}[ex]
def call(i):
    b = n.NCycDB_best[i]
    if not t.Contig[i]:
        return 'n.a.'
    if not b:
        return 'no hit'
    if b.startswith('(homolog'):
        return 'homolog outside the gene map'
    ok = float(n.NCycDB_identity[i] or 0) >= 40 and float(n.NCycDB_ref_cov[i] or 0) >= 70
    return b if ok else f'{b} (below threshold)'
def where(i):
    b, c = t.Bin[i], t.Contig[i]
    return f'{b} ({c})' if c else b
S2B = pd.DataFrame({'Gene': t.Gene, 'Bin or contig': [where(i) for i in t.index], 'Reference accession': t['Reference (accession)'],
                    'Amino acid identity (%)': t['Identity (%)'].map(r1), 'Coverage of the reference (%)': t['Reference coverage (%)'].map(r1),
                    'Evidence level': t.Evidence, 'DRAM/KOfam detection': [dram(t.loc[i]) for i in t.index],
                    'NCycDB call (gene family)': [call(i) for i in t.index], 'NCycDB status': [status(i) for i in t.index], 'Note': t.Note})
sup = S2B[S2B['Evidence level'] == 'supported']
counts = sup['NCycDB status'].value_counts()
assert len(sup) == 112 and counts['out of scope'] == 19 and counts['concordant'] == 47 and (sup['NCycDB status'] != 'n.a.').all()

# ------------------------------------------------------------------ S2D: Anammox panel best hits (seção 29)
a = pd.read_csv(f'{RC}/29_panel_masking0/orf_assignment_masking0.tsv', sep='\t')
bl = pd.read_csv(f'{RC}/06_tableS2/hzs_orfs_vs_anammox_panel_long.tsv', sep='\t', header=None,
                 names=['orf', 'acc', 'pident', 'alnlen', 'qs', 'qe', 'ss', 'se', 'evalue', 'bitscore', 'qlen', 'slen'])
bl = bl.sort_values('bitscore', ascending=False).drop_duplicates(['orf', 'acc']).set_index(['orf', 'acc']).pident
porg = {x.accession: HDR.get(x.accession, '') for _, x in panel.iterrows()}
S2D = pd.DataFrame({'Contig': a.contig, 'ORF': a.orf, 'Length (aa)': a.orf_len_aa, 'Partial (yes/no)': a.partial,
                    'Subunit assigned by sequence': a.subunit_by_sequence,
                    '5C2V chain identity (%)': [f'{"beta" if s == "hzsB" else "gamma"} {r1(p)}' if pd.notna(p) else 'n.a.'
                                                for s, p in zip(a.chain_5C2V, a.chain_5C2V_pident)],
                    'Best Brocadia hit (accession)': a.brocadia_ref, 'Best Brocadia hit (organism)': a.brocadia_ref.map(porg),
                    'GenBank label corrected': a.brocadia_ref.map(lambda x: f'yes (GenBank label {"hzsC" if FIX[x] == "hzsB" else "hzsB"})' if x in FIX else 'no'),
                    'Amino acid identity (%) (--masking 0)': a.brocadia_pident.map(r1),
                    'Best Kuenenia hit (accession)': a.kuenenia_ref, 'Best Kuenenia hit identity (%)': a.kuenenia_pident.map(r1),
                    'Margin (points)': a.margin_vs_kuenenia.map(r1), 'Used in range': (a.orf_len_aa >= 100).map({True: 'yes', False: 'no'}),
                    'BLASTp identity (%) (06_tableS2)': [r1(bl.get((o, r))) if (o, r) in bl.index else 'n.a. (reference not in the BLASTp set)'
                                                         for o, r in zip(a.orf, a.brocadia_ref)]})
GO = {'hzsA': 0, 'hzsB': 1, 'hzsC': 2, 'hdh': 3}
S2D = S2D.sort_values(['Subunit assigned by sequence', 'Contig', 'ORF'], key=lambda s: s.map(GO) if s.name == 'Subunit assigned by sequence' else s).reset_index(drop=True)

# ------------------------------------------------------------------ S2C: Notes
cat = ['concordant', 'same gene below threshold', 'pmo (amo-assigned)', 'nar/nxr family', 'hzsB/hzsC nomenclature swap',
       'same-family homolog outside the map', 'other gene below threshold', 'out of scope']
N = [('Sheet S2A', 'References', 'One row per reference sequence (accession). Gene/protein: gene searched or represented. NCBI record organism: organism as given in the NCBI/UniProt record (RefSeq MULTISPECIES records give a genus or higher rank). Reference organism cited: organism used for this reference in the text (spelling normalized; equal to the record organism when the reference is not named in the text). Use: targeted search (BLASTp against bin proteomes and tBLASTn against the assembly); anammox panel (22 sequences used to assign the genus of the hzs/hdh contigs; sheet S2D); HzsA tree (12 references of the HzsA phylogeny); DsrA oxidative / DsrA reductive (references of the DsrA typing tree); AmoA tree (Figure S3); other. Source database and length (aa) as in the record. Accession versions are given as in the files used.'),
     ('Sheet S2B', 'Gene; Bin or contig', 'One row per gene and bin or contig (392 rows). Bin or contig: consolidated bin and contig carrying the best targeted hit; "unbinned" = contig not in any of the 25 bins; ‡ = unbinned contig attributed to the bin by depth, GC, tetranucleotide frequency, co-binning and read pairs; a bin without a contig = gene not detected in that bin.'),
     ('Sheet S2B', 'Reference accession; Amino acid identity (%); Coverage of the reference (%)', 'Reference with the highest bit score at the locus (sheet S2A). Identity: BLASTp on the predicted open reading frame for contigs of the 25 bins, tBLASTn against the assembly for unbinned contigs. Coverage of the reference: aligned fraction of the reference, summed over adjacent open reading frames when a gene was split. Values with one decimal.'),
     ('Sheet S2B', 'Evidence level', 'supported: E-value <=1e-10, amino acid identity >=40% and coverage >=70% of the reference. partial: same identity and E-value criteria with coverage <70% of the reference. KOfam only: gene annotated by DRAM/KOfam in a bin without a targeted hit >=40% identity at the locus. not detected: no targeted hit >=40% identity and no KOfam call in the bin. Close paralogs (hdh and hao, napA and formate dehydrogenase, narG and nxrA, amo and pmo) were separated by DRAM and Pfam annotation, reciprocal best hit, or host-bin taxonomy and coverage (Supplementary Methods S1).'),
     ('Sheet S2B', 'DRAM/KOfam detection', 'yes: the DRAM gene at the locus carries the expected KO; partial: domain-level annotation without KO assignment (PF00384 for narG/nxrA); no: no corresponding annotation, or contig outside the consolidated bins (unbinned contigs have no DRAM annotation).'),
     ('Sheet S2B', 'NCycDB call (gene family)', 'Gene family of the best hit of the open reading frame at the locus in the NCycDB 95% identity representative set (DIAMOND blastp, --more-sensitive, E <=1e-10; 206,931 open reading frames of the assembly). A call requires identity >=40% and coverage >=70% of the NCycDB reference; otherwise the family is given with "(below threshold)". "homolog outside the gene map": best hit to a sequence not mapped to a gene family; "no hit": no hit at E <=1e-10; n.a.: no locus (gene not detected).'),
     ('Sheet S2B', 'NCycDB status', 'Assigned only to the 112 supported rows (n.a. in the other rows). ' + ' '.join([
         'concordant: NCycDB calls the same gene.', 'same gene below threshold: the best NCycDB hit is the same gene but below the identity or coverage threshold.',
         'pmo (amo-assigned): NCycDB calls pmoA/B/C for a gene assigned to amoA/B/C; amo and pmo belong to the same copper membrane monooxygenase family.',
         'nar/nxr family: NCycDB calls another gene of the NarG/NxrA family (narG, narZ, narH/nxrB); NCycDB groups nitrite oxidoreductase with nitrate reductase.',
         'hzsB/hzsC nomenclature swap: NCycDB (as TrEMBL) uses beta and gamma swapped relative to Swiss-Prot and PDB 5C2V, which this table follows.',
         'same-family homolog outside the map: best hit to a homolog of the same family that is not in the NCycDB gene map.',
         'other gene below threshold: best hit to another gene of the same family (nar/nxr or pmo) below the threshold.',
         'out of scope: sox, dsr and mcr genes, which NCycDB does not cover.'])),
     ('Sheet S2D', 'Anammox panel best hits', 'One row per open reading frame (Prodigal v2.6.3, meta mode) with an hzs or hdh hit on the ten unbinned Ca. Brocadia-affiliated contigs, aligned with DIAMOND v2.1.24 (blastp, --more-sensitive, E <=1e-10, --masking 0, --comp-based-stats 0) to the 22-sequence anammox panel (sheet S2A, Use = anammox panel). Subunit assigned by sequence: subunit of the panel reference with the highest bit score, with GenBank labels corrected, checked against the beta and gamma chains of PDB 5C2V (Q1Q0T4, Q1Q0T3; BLASTp) and against gene order in the operon (hzsB-hzsC-hzsA in all four contigs with more than one hzs gene). Best Brocadia hit and best Kuenenia hit: highest identity among the panel references of that genus for the same subunit. Margin: Brocadia minus Kuenenia identity (points). Used in range: yes for open reading frames >=100 aa; the two 56-aa fragments (NODE_5668_4 and NODE_17695_1) are listed but not used. Partial: open reading frame truncated at a contig end (Prodigal). BLASTp identity: BLASTp of the same open reading frame against the same reference (n.a. when the reference was retrieved from NCBI for the panel and not included in the BLASTp set). NODE_24575 carries only a 47-aa tBLASTn fragment of hzsB (sheet S2B) and no open reading frame with a panel hit.'),
     ('Sheet S2D', 'GenBank label swap', 'The GenBank records of Ca. Brocadia fulgida carry swapped labels: KKO20885.1 ("hydrazine synthase subunit C") is HZS-beta (hzsB) and KKO20886.1 ("hydrazine synthase subunit B") is HZS-gamma (hzsC) by sequence identity to the Swiss-Prot entries Q1Q0T4 (HZS-beta) and Q1Q0T3 (HZS-gamma), the beta and gamma chains of the hydrazine synthase structure PDB 5C2V. This table uses the corrected subunits.'),
     ('Sheet S2D', '--masking 0', 'DIAMOND v2.1.24 masks low-complexity regions (tantan) and applies composition-based statistics by default; with these defaults the identity of NODE_7170_3 (hzsB, a beta-propeller with internal repeats) was 10-18 points lower than BLASTp over the same aligned positions. Masking and composition-based statistics were therefore turned off (--masking 0, --comp-based-stats 0), which reproduces the BLASTp identities of the best hits.'),
     ('All sheets', 'Targeted-search criterion', 'A gene was considered supported when the best hit to the curated reference had an E-value <=1e-10, an amino acid identity >=40%, and a coverage of at least 70% of the reference sequence; hits covering less were reported as partial. Identities and coverages are given with one decimal.')]
S2C = pd.DataFrame(N, columns=['Sheet', 'Item', 'Description'])

# ------------------------------------------------------------------ escrita
sheets = {'S2A References': S2A, 'S2B Targeted search NCycDB': S2B, 'S2C Notes': S2C, 'S2D Anammox panel best hits': S2D}
csvname = {'S2A References': 'Table_S2A_References.csv', 'S2B Targeted search NCycDB': 'Table_S2B_Targeted_search_NCycDB.csv',
           'S2C Notes': 'Table_S2C_Notes.csv', 'S2D Anammox panel best hits': 'Table_S2D_Anammox_panel_best_hits.csv'}
for nm, df in sheets.items():
    for c in df.columns:
        assert not df[c].astype(str).str.startswith('=').any(), (nm, c)   # openpyxl grava '=...' como fórmula
    df.to_csv(f'{OUT}/{csvname[nm]}', index=False, encoding='utf-8')
xl = f'{OUT}/Table_S2.xlsx'
with pd.ExcelWriter(xl, engine='openpyxl') as w:
    for nm, df in sheets.items():
        df.to_excel(w, sheet_name=nm, index=False)
wb = load_workbook(xl)
for ws in wb.worksheets:
    for c in ws[1]:
        c.font = Font(bold=True)
    for col in ws.columns:
        L = max(len(str(c.value or '')) for c in col)
        ws.column_dimensions[col[0].column_letter].width = min(max(10, L + 2), 60 if ws.title != 'S2C Notes' else 110)
        num = '(%)' in str(col[0].value) or str(col[0].value).startswith(('Margin', 'Length'))
        for c in col[1:]:
            if num and isinstance(c.value, str) and re.fullmatch(r'-?\d+(\.\d)?', c.value):
                c.value = float(c.value); c.number_format = '0' if str(col[0].value).startswith('Length') else '0.0'
            elif num and isinstance(c.value, (int, float)):
                c.number_format = '0' if str(col[0].value).startswith('Length') else '0.0'
            c.alignment = Alignment(wrap_text=ws.title == 'S2C Notes', vertical='top')
    ws.freeze_panes = 'A2'
p = wb.properties
p.creator = ''; p.lastModifiedBy = ''; p.title = 'Table S2'; p.company = None if hasattr(p, 'company') else None
wb.save(xl)

print({nm: len(df) for nm, df in sheets.items()})
print(counts.reindex(cat).to_string())
print('S2A Use:', S2A.Use.str.split('; ').explode().value_counts().to_dict())
