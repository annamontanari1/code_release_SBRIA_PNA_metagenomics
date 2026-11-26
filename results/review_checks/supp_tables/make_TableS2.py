"""Table S2 (SBRIA shotgun): um acerto de busca dirigida por linha + linhas KOfam only / not detected; aba References.

Critério (Methods): e-value <= 1e-10; identidade >= 40%; cobertura da referência >= 70% -> supported; < 70% -> partial.
Famílias com parálogos (hdh/hao, napA/formato desidrogenase, narG/nxrA): o acerto só é aceito se a anotação DRAM/Pfam do gene
sobreposto for concordante; sem gene DRAM no locus (contig não binado ou < 2.500 pb) o acerto é aceito com nota.
narG x nxrA: função atribuída pela taxonomia do bin hospedeiro (Nitrospira e Brocadia -> nxrA; demais -> narG).
Cada locus (ORF, ou HSP quando não há ORF) é atribuído ao gene de maior bit score; uma linha por gene x contig (melhor acerto).
Não altera arquivos existentes: lê results/, DRAM annotation/, das_tool_result (1)/ e results/review_checks/; escreve em supp_tables/.
"""
import os, re, csv, collections, itertools, random
import numpy as np, pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment

HERE = os.path.dirname(os.path.abspath(__file__))
RC = os.path.dirname(HERE); ROOT = os.path.abspath(os.path.join(RC, '..', '..'))
D = os.path.join(ROOT, 'das_tool_result (1)')
EV, ID, COV = 1e-10, 40.0, 70.0
node = lambda x: re.match(r'(NODE_\d+)', x).group(1)

# ---------------- contigs: comprimento, GC, profundidade, bins ----------------
INFO = {}
for l in open(os.path.join(RC, '00_base', 'contig_len_gc.tsv')):
    n, L, g = l.rstrip('\n').split('\t'); INFO[node(n)] = dict(len=int(L), gc=float(g))
for r in csv.reader(open(os.path.join(ROOT, 'results', 'contig_check', 'asm_coverage.tsv')), delimiter='\t'):
    if not r[0].startswith('#'): INFO[node(r[0])]['depth'] = float(r[6])
KEY = {r[0].replace('.fa', ''): r[1].replace('.fasta', '') for r in csv.reader(open(os.path.join(D, 'das_tool_output_dir_DASTool_bins', 'das_tool_name_key.tsv')), delimiter='\t') if not r[0].startswith('Original')}
FINAL = {node(l.split()[0]): KEY[l.split()[1]] for l in open(os.path.join(D, 'das_tool_output_dir_DASTool_scaffolds2bin.txt'))}
BINNER = collections.defaultdict(dict)
for t, nm in [('metabat2', 'MetaBAT2'), ('maxbin2', 'MaxBin2'), ('concoct', 'CONCOCT')]:
    for l in open(os.path.join(D, f'NGS1140_6S_{t}_bins_contigs_to_bins.tsv')):
        c, b = l.split()[:2]; BINNER[node(c)][nm] = b.replace(f'NGS1140_6S_{t}_bins.', '')
def binners(n):
    d = BINNER.get(n, {}); return '; '.join(f'{k} {d[k]}' for k in ['MetaBAT2', 'MaxBin2', 'CONCOCT'] if k in d) or 'none'
TAX = {}
for f in ['bacteria.csv', 'archea.csv']:
    for r in csv.DictReader(open(os.path.join(ROOT, f))): TAX[r['User Genome'].replace('.fa', '')] = r['Classification']
def genus(b):
    t = TAX.get(b, '')
    for r in ['g__', 'f__', 'o__']:
        for x in t.split(';'):
            if x.startswith(r) and x[3:]: return x[3:]
    return ''

# ---------------- DRAM: genes por contig (coordenadas, KO, Pfam) ----------------
DRAM = collections.defaultdict(list)
for r in csv.reader(open(os.path.join(ROOT, 'DRAM annotation', 'annotations.tsv')), delimiter='\t'):
    if not r[0] or r[1] == 'fasta': continue
    DRAM[node(r[2])].append(dict(id=node(r[2]) + '_' + r[0].rsplit('_', 1)[1], bin=r[1].split('.fasta')[0], s=int(r[4]), e=int(r[5]),
                                 strand='+' if r[6] == '1' else '-', ko=r[8], desc=r[9], pfam=r[17]))
def dram_overlap(n, s, e):
    s, e = min(s, e), max(s, e)
    ov = [(min(e, g['e']) - max(s, g['s']), g) for g in DRAM.get(n, []) if min(e, g['e']) > max(s, g['s'])]
    return max(ov, key=lambda x: x[0])[1] if ov else None

# ---------------- referências e genes ----------------
SEQLEN, DESC = {}, {}
def readfa(fn):
    k = None
    for l in open(fn):
        if l[0] == '>':
            a = l[1:].strip().split(' ', 1)
            ids = [x for x in a[0].split('|') if re.match(r'^[A-Z0-9_]+\.\d+$|^Q1PZD8$|^8SR2_C$', x)]
            k = ids[0] if ids else a[0]
            # não sobrescreve: amo_refs.fasta (AmoB/AmoC truncadas em 70 aa) é lida depois das versões oficiais
            k = k if k not in SEQLEN else None
            if k: SEQLEN[k] = 0; DESC.setdefault(k, a[1] if len(a) > 1 else '')
        elif k: SEQLEN[k] += len(l.strip())
TS = os.path.join(ROOT, 'results', 'targeted_search')
for f in [os.path.join(RC, 'refs', 'refs_clean.faa'), os.path.join(HERE, 'work', 'extra_refs.faa'), os.path.join(RC, 'refs', 'pmoC_refs.faa'),
          os.path.join(RC, 'refs', 'dsrA_refs_panel.faa'), os.path.join(TS, 'pseudomonas_nirS.fasta'), os.path.join(TS, 'scalindua_nirS.fasta'),
          os.path.join(TS, 'hao_ref.fasta'), os.path.join(TS, 'amo_refs.fasta'), os.path.join(TS, 'nxr_anammox_ref.fasta')]:
    if os.path.exists(f): readfa(f)

GENES = {'amoA': ['Q04507.2', 'WP_062484767.1', 'TKS59187.1'], 'amoB': ['Q04508.2', 'WP_062484768.1', 'TKS59188.1'],
         'amoC': ['WP_011111554.1', 'WP_062484140.1', 'TKS59186.1'],
         'narG/nxrA': ['YEH58028.1', 'Q1PZD8', 'WP_230999932.1'], 'nxrB': ['YEH58030.1'],
         'hzsA': ['WP_099325812.1'], 'hzsB': ['WP_164994616.1'], 'hzsC': ['WP_478210066.1'], 'hdh': ['Q1PW30.2'],
         'napA': ['WP_104493079.1'], 'nirS': ['WP_164995172.1', 'Q51700.1'], 'nirK': ['WP_316989695.1'], 'nrfA': ['P0ABK9.1'],
         'nirB': ['P08201.4'], 'nirD': ['P0A9I8.1'], 'norB': ['BAA32546.1'], 'norC': ['Q51662.3'],
         # nosZ: clado I (Paracoccus) e clado II (Ca. D. denitrificans/Okubo, Thiothrix, Dechloromonas, Wolinella, Anaeromyxobacter);
         # o melhor bit score define a referência de cada acerto
         'nosZ': ['WP_230999923.1', 'BBO22106.1', 'WP_301017512.1', 'AAZ46324.1', 'WP_129545366.1', 'WP_011421454.1'], 'mcrA': ['WP_505736924.1'],
         'soxA': ['WP_024845877.1'], 'soxB': ['WP_024845878.1'], 'soxC': ['WP_024845879.1'], 'soxD': ['WP_024845880.1'],
         'soxX': ['WP_024845874.1'], 'soxY': ['WP_024845875.1'], 'soxZ': ['WP_024845876.1'],
         'dsrA': ['BBO20425.1'], 'dsrB': ['BBO20424.1']}
REF2GENE = {r: g for g, rs in GENES.items() for r in rs}
KO = {'amoA': 'K10944', 'amoB': 'K10945', 'amoC': 'K10946', 'narG/nxrA': 'K00370', 'nxrB': 'K00371', 'hdh': 'K20935', 'napA': 'K02567',
      'nirS': 'K15864', 'nirK': 'K00368', 'nrfA': 'K03385', 'nirB': 'K00362', 'nirD': 'K00363', 'norB': 'K04561', 'norC': 'K02305',
      'nosZ': 'K00376', 'mcrA': 'K00399', 'soxA': 'K17222', 'soxB': 'K17224', 'soxC': 'K17225', 'soxD': 'K22622', 'soxX': 'K17223',
      'soxY': 'K17226', 'soxZ': 'K17227', 'dsrA': 'K11180', 'dsrB': 'K11181'}
PARALOG_FAMILIES = {'hdh': ({'K20935'}, set()), 'napA': ({'K02567'}, set()), 'narG/nxrA': ({'K00370'}, {'PF00384'})}

# ---------------- acertos tBLASTn (contig) e BLASTp (ORF prodigal) ----------------
def load(fn, kind):
    out = []
    for l in open(fn):
        p = l.rstrip('\n').split('\t')
        if float(p[10]) > EV or p[0] not in REF2GENE: continue
        out.append(dict(ref=p[0], subj=p[1], contig=node(p[1]), pid=float(p[2]), aln=int(p[3]), qs=int(p[6]), qe=int(p[7]),
                        ss=int(p[8]), se=int(p[9]), ev=float(p[10]), bits=float(p[11]), qlen=int(p[12]), slen=int(p[13]), kind=kind))
    return out
TBN = load(os.path.join(RC, '00_base', 'tblastn_allrefs_vs_assembly.tsv'), 'tblastn') + load(os.path.join(HERE, 'work', 'tblastn_sox_dsr_vs_assembly.tsv'), 'tblastn')
BP = load(os.path.join(RC, '00_base', 'blastp_allrefs_vs_hitORFs.tsv'), 'blastp') + load(os.path.join(HERE, 'work', 'blastp_sox_dsr_vs_ORFs.tsv'), 'blastp')
ORF = {}
for fn in [os.path.join(RC, '00_base', 'hit_contigs_orfs.faa'), os.path.join(HERE, 'work', 'c_orfs.faa')]:
    for l in open(fn):
        if l[0] == '>':
            p = l[1:].split(' # '); ORF[p[0]] = (int(p[1]), int(p[2]), '+' if p[3] == '1' else '-', re.search(r'partial=(\d\d)', p[4]).group(1))
def union(iv):
    tot, cs, ce = 0, None, None
    for s, e in sorted(iv):
        if cs is None or s > ce + 1:
            if cs is not None: tot += ce - cs + 1
            cs, ce = s, e
        else: ce = max(ce, e)
    return tot + (ce - cs + 1 if cs is not None else 0)

cand = {}   # (gene, contig) -> melhor referência e métricas
for g, refs in GENES.items():
    for ref in refs:
        for c in sorted(set(h['contig'] for h in TBN if h['ref'] == ref)):
            hs = [h for h in TBN if h['ref'] == ref and h['contig'] == c]
            best = max(hs, key=lambda h: h['bits'])
            ucov = 100 * union([(h['qs'], h['qe']) for h in hs]) / SEQLEN[ref]
            bps = [h for h in BP if h['ref'] == ref and h['contig'] == c]
            ob = max(bps, key=lambda h: h['bits']) if bps else None
            bits = max(best['bits'], ob['bits'] if ob else 0)
            k = (g, c)
            if k in cand and cand[k]['bits'] >= bits: continue
            cand[k] = dict(gene=g, contig=c, ref=ref, tbn=best, ucov=ucov, orf=ob, bits=bits,
                           others={r: max([h['pid'] for h in TBN if h['ref'] == r and h['contig'] == c] or [None], key=lambda x: x or 0) for r in refs if r != ref})

# locus -> gene de maior bit score (evita contar o mesmo ORF em dois genes)
locus_owner = {}
for k, v in cand.items():
    loc = v['orf']['subj'] if v['orf'] else f"{v['contig']}:{min(v['tbn']['ss'], v['tbn']['se']) // 500}"
    v['locus'] = loc
    if loc not in locus_owner or cand[locus_owner[loc]]['bits'] < v['bits']: locus_owner[loc] = k

# ---------------- TNF/GC contra bin.009 (hzs/hdh) ----------------
SEQ = {}
need = {c for (g, c) in cand if g in ('hzsA', 'hzsB', 'hzsC', 'hdh')} | {c for c, b in FINAL.items() if b == 'bin.009'}
k = None
for l in open(os.path.join(ROOT, 'results', 'targeted_search', 'bin_assemblies', 'NGS1140_6S_assembly.fa')):
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
B9 = [c for c, b in FINAL.items() if b == 'bin.009']
C9 = sum(tnf(SEQ[c]) * len(SEQ[c]) for c in B9); C9 /= C9.sum()
b9s = [SEQ[c] for c in B9]; b9w = [len(s) for s in b9s]
GC9 = np.array([INFO[c]['gc'] for c in B9]); GCM, GCSD = GC9.mean(), GC9.std(ddof=1)
def tnf_test(c):
    rng = random.Random(int(c.split('_')[1]))   # semente por contig: resultado independe da ordem de execução
    s = SEQ[c]; d = np.linalg.norm(tnf(s) - C9); L = len(s); nd = []
    while len(nd) < 500:
        x = rng.choices(b9s, weights=b9w)[0]
        if len(x) >= L: i = rng.randrange(len(x) - L + 1); nd.append(np.linalg.norm(tnf(x[i:i + L]) - C9))
    nd = np.array(nd); return d, np.percentile(nd, 95), 100 * (nd < d).mean()

# ---------------- montar linhas ----------------
NITRO = {'bin.004', 'bin.007'}
# amo -> pmo pelo NCycDB (14_ncycdb/TableS2_vs_NCycDB.tsv; followup_21out, 21/10/2026). PMO_NCYC: melhor hit pmoA/B/C com
# 84.4-98.8% de identidade e cobertura da referência >= 70%; PMO_NCYC_LOWCOV: melhor hit pmo com cobertura < 70%.
PMO_NCYC = {('amoA', 'NODE_18249'), ('amoA', 'NODE_7484'), ('amoA', 'NODE_22767'), ('amoA', 'NODE_11718'),
            ('amoB', 'NODE_18249'), ('amoB', 'NODE_11718'),
            ('amoC', 'NODE_17726'), ('amoC', 'NODE_2742'), ('amoC', 'NODE_7852'), ('amoC', 'NODE_12015'), ('amoC', 'NODE_10604')}
PMO_NCYC_LOWCOV = {('amoA', 'NODE_7852'), ('amoB', 'NODE_22767'), ('amoB', 'NODE_7484'),
                   ('amoC', 'NODE_32378'), ('amoC', 'NODE_17747'), ('amoC', 'NODE_18653'), ('amoC', 'NODE_11718')}
# nxrB: notas do painel NxrB/NarH (14_ncycdb/followup_21out/nxrB_orfs_vs_panel.tsv); NODE_835/NODE_136 -> "nxrB-like"
NXRB_NOTE = {'NODE_16231': "Brocadia-type NxrB (96.5% B. sinica, 96.0% B. fulgida); NCycDB 'narH' family contains anammox NxrB",
             'NODE_5446': 'NxrB, Nitrospira/anammox clade (63%); host unresolved',
             'NODE_835': 'beta subunit of periplasmic NXR clade in a non-nitrifying host; function (NO2- oxidation vs NO3- reduction) not resolved by sequence',
             'NODE_136': 'beta subunit of periplasmic NXR clade in a non-nitrifying host; function (NO2- oxidation vs NO3- reduction) not resolved by sequence'}
NXRB_LIKE = {'NODE_835', 'NODE_136'}
def fmt_e(x): return f'{x:.1e}' if x > 0 else '0.0'
def r1(x): return round(float(x) + 1e-9, 1)
rows, excluded = [], []
for k, v in cand.items():
    if locus_owner[v['locus']] != k: continue
    g, c, t, o = v['gene'], v['contig'], v['tbn'], v['orf']
    # Methods: BLASTp contra os proteomas dos 25 bins; tBLASTn contra o assembly completo.
    # Identidade principal: BLASTp no ORF para contigs dos bins; tBLASTn para contigs não binados. As duas ficam em colunas próprias.
    orf_pid = o['pid'] if o else None
    tbn_pid = t['pid']
    pid = orf_pid if (c in FINAL and o) else tbn_pid
    ev = min(t['ev'], o['ev']) if o else t['ev']
    refcov = min(100.0, max(v['ucov'], 100 * (o['qe'] - o['qs'] + 1) / SEQLEN[v['ref']] if o else 0))
    qcov = 100 * (o['se'] - o['ss'] + 1) / o['slen'] if o else None
    if o: s, e, strand, part = ORF[o['subj']]
    else: s, e = sorted([t['ss'], t['se']]); strand = '+' if t['se'] > t['ss'] else '-'; part = ''
    b = FINAL.get(c, 'unbinned'); I = INFO[c]
    gname = g
    if g == 'narG/nxrA':
        gname = 'nxrA' if (b in NITRO | {'bin.009'} or c == 'NODE_2432' or (b == 'unbinned' and v['ref'] in ('YEH58028.1', 'Q1PZD8') and pid >= 70)) else ('narG' if b != 'unbinned' else 'narG/nxrA (unassigned)')
    dg = dram_overlap(c, s, e)
    notes = []
    status = None
    if pid < ID:
        if dg and dg['ko'] == KO.get(g):
            status = 'KOfam only'; notes.append(f'best targeted hit {pid:.1f}% identity (<40%)')
        else:
            continue
    if status is None and g in PARALOG_FAMILIES:
        conc, pf = PARALOG_FAMILIES[g]
        if dg is None: notes.append('no DRAM gene at locus (unbinned or contig <2,500 bp); paralog check not possible')
        elif dg['ko'] and dg['ko'] not in conc:
            excluded.append(dict(gene=gname, contig=c, bin=b, ref=v['ref'], identity=round(pid, 1), ref_cov=round(refcov, 1), DRAM_KO=dg['ko'], DRAM_desc=dg['desc'][:60], reason=f"DRAM KO at locus discordant ({dg['ko']})")); continue
        elif not dg['ko'] and not any(p in dg['pfam'] for p in pf):
            notes.append('DRAM gene without KO/concordant Pfam')
        else:
            notes.append(f"annotation concordant ({dg['ko'] or 'Pfam ' + ','.join(p for p in pf if p in dg['pfam'])})")
    if status is None:
        status = 'supported' if refcov >= COV else 'partial'
    if part in ('10', '01', '11'): notes.append('ORF truncated at contig end')
    other = [f"{r} {p:.1f}%" for r, p in v['others'].items() if p]
    if other: notes.append('other refs: ' + '; '.join(other))
    if g in ('amoA', 'amoB', 'amoC'):
        if b == 'bin.014': notes.append('pmo of Methylocystis (bin.014), not amo')
        elif c == 'NODE_4969': notes.append('Nitrosomonas amoCAB operon (AmoA 99.3% N. oligotropha)')
        elif c == 'NODE_13511': notes.append('comammox clade A amoAB (92.9% Ca. N. inopinata AmoA)')
        elif c == 'NODE_2578': notes.append('comammox-type amoC (85.8% Ca. N. inopinata)')
        elif (g, c) in PMO_NCYC: notes.append('pmo, not amo (NCycDB pmoA/B/C 84-99%)')
        elif (g, c) in PMO_NCYC_LOWCOV: notes.append('closer to pmo than to amo (NCycDB, ref coverage < 70%); affiliation unresolved')
        elif b == 'unbinned': notes.append('amo/pmo affiliation unresolved')
    if dg and dg['bin'] != b and b != 'unbinned': notes.append(f"DRAM bin {dg['bin']}")
    row = dict(gene=gname, bin=b, bin_taxon=genus(b) if b != 'unbinned' else '', contig=c, contig_length_bp=I['len'], depth_x=r1(I['depth']),
               contig_GC_pct=r1(I['gc']), hit_coordinates=f'{s}-{e} ({strand})', reference=v['ref'], identity_pct=r1(pid), tbn_identity_pct=r1(tbn_pid), orf_identity_pct=r1(orf_pid) if orf_pid is not None else '',
               query_coverage_pct=r1(qcov) if qcov is not None else '', reference_coverage_pct=r1(refcov), evalue=fmt_e(ev), bit_score=r1(v['bits']),
               evidence=status, binner=binners(c), DAS_Tool_retained='yes' if c in FINAL else 'no', DRAM_gene=dg['id'] if dg else '',
               DRAM_KO=dg['ko'] if dg else '', note='; '.join(notes))
    if g in ('hzsA', 'hzsB', 'hzsC', 'hdh'):
        d, p95, pct = tnf_test(c)
        row.update(GC_z_vs_bin009=round((I['gc'] - GCM) / GCSD, 2), TNF_dist_to_bin009=round(d, 4), TNF_null_P95_bin009=round(p95, 4), TNF_percentile_vs_bin009=r1(pct))
    rows.append(row)

# KOfam only: genes-alvo anotados pelo DRAM sem acerto dirigido aceito no locus
acc = collections.defaultdict(list)
for r in rows:
    s, e = map(int, r['hit_coordinates'].split(' ')[0].split('-')); acc[r['contig']].append((s, e))
INV = {v: k for k, v in KO.items()}
for c, gl in DRAM.items():
    for gd in gl:
        if gd['ko'] not in INV: continue
        if any(min(e, gd['e']) > max(s, gd['s']) for s, e in acc.get(c, [])): continue
        g = INV[gd['ko']]; b = gd['bin']; I = INFO[c]
        gname = ('nxrA' if b in NITRO | {'bin.009'} else 'narG') if g == 'narG/nxrA' else g
        best = [h for h in TBN if REF2GENE[h['ref']] == g and h['contig'] == c and min(max(h['ss'], h['se']), gd['e']) > max(min(h['ss'], h['se']), gd['s'])]
        bh = max(best, key=lambda h: h['bits']) if best else None
        same = [r for r in rows if r['contig'] == c and r['gene'] == gname and r['evidence'] in ('supported', 'partial')]
        if bh and bh['pid'] >= ID:
            # cópia adicional / ORF dividido no mesmo contig: classificado pelo critério, não como KOfam only
            rc = min(100.0, 100 * union([(h['qs'], h['qe']) for h in best if h['ref'] == bh['ref']]) / SEQLEN[bh['ref']])
            ev_ = 'supported' if rc >= COV else 'partial'
            nt = ('additional copy or split ORF on a contig with another accepted hit' if same else 'locus not captured by the per-contig best hit')
            nt += '; annotation concordant (' + gd['ko'] + ')'
            rows.append(dict(gene=gname, bin=b, bin_taxon=genus(b), contig=c, contig_length_bp=I['len'], depth_x=r1(I['depth']), contig_GC_pct=r1(I['gc']),
                             hit_coordinates=f"{gd['s']}-{gd['e']} ({gd['strand']})", reference=bh['ref'], identity_pct=r1(bh['pid']), query_coverage_pct='',
                             reference_coverage_pct=r1(rc), evalue=fmt_e(bh['ev']), bit_score=r1(bh['bits']), evidence=ev_, binner=binners(c),
                             DAS_Tool_retained='yes', DRAM_gene=gd['id'], DRAM_KO=gd['ko'], note=nt))
            acc[c].append((gd['s'], gd['e'])); continue
        note = (f"best targeted hit {bh['pid']:.1f}% ({bh['ref']}; {100 * bh['aln'] / SEQLEN[bh['ref']]:.0f}% of reference)" if bh else 'no targeted hit (e<=1e-10)')
        if g in PARALOG_FAMILIES: note += '; paralog family: KOfam call not confirmed'
        rows.append(dict(gene=gname, bin=b, bin_taxon=genus(b), contig=c, contig_length_bp=I['len'], depth_x=r1(I['depth']), contig_GC_pct=r1(I['gc']),
                         hit_coordinates=f"{gd['s']}-{gd['e']} ({gd['strand']})", reference=bh['ref'] if bh else '', identity_pct=r1(bh['pid']) if bh else '',
                         query_coverage_pct='', reference_coverage_pct=r1(100 * bh['aln'] / SEQLEN[bh['ref']]) if bh else '', evalue=fmt_e(bh['ev']) if bh else '',
                         bit_score=r1(bh['bits']) if bh else '', evidence='KOfam only', binner=binners(c), DAS_Tool_retained='yes',
                         DRAM_gene=gd['id'], DRAM_KO=gd['ko'], note=note))

# Regra HAO (19/10/2026): contig com acerto "hdh" cujo melhor hit recíproco no proteoma de K. stuttgartiensis CSTR1
# (NZ_CP049055.1) é HAO (KsCSTR_RS17015) ou o citocromo multi-heme da família HAO (KsCSTR_RS16015), com diferença
# >= 42 pontos sobre a identidade às HDH (KsCSTR_RS18435/RS04630) -> parálogo, movido para a aba de excluídos.
HAO = {}
for l in open(os.path.join(RC, 'comammox', 'hdh_reciprocal', 'hdh_vs_HAO_family_assignment.tsv')):
    if l.startswith(('#', 'contig')): continue
    p = l.rstrip('\n').split('\t')
    if p[-1] == 'HAO-family paralog' and float(p[8]) >= 42: HAO[p[0]] = (p[5], p[6], p[7], p[8])
keep = []
for r in rows:
    if r['gene'] == 'hdh' and r['contig'] in HAO:
        fam, pb, ph, mg = HAO[r['contig']]
        excluded.append(dict(gene='hdh', contig=r['contig'], bin=r['bin'], ref=r['reference'], identity=r['identity_pct'],
                             ref_cov=r['reference_coverage_pct'], DRAM_KO=r['DRAM_KO'], DRAM_desc='',
                             reason=f"best reciprocal hit in K. stuttgartiensis NZ_CP049055.1 is HAO-family ({fam}, {float(pb):.1f}%) vs HDH {float(ph):.1f}% (margin {float(mg):.1f} points >= 42)"))
    else: keep.append(r)
rows = keep
# NODE_2578 (amoC, bin.007) e NODE_13511 (amoA/amoB, não binado) atribuídos ao bin.007 (verificação comammox, 19/10/2026)
IND = 'attributed to bin.007 by depth, GC, TNF, MaxBin2 co-binning and 46 read pairs'
for r in rows:
    if r['contig'] == 'NODE_13511' and r['gene'] in ('amoA', 'amoB'):
        r['bin'] = 'bin.007‡'; r['bin_taxon'] = 'Nitrospira_F'; r['evidence'] = 'supported'
        r['note'] = IND + ' (‡ indirect attribution; contig not retained by DAS Tool); comammox clade A (AmoA 92.9% Ca. N. inopinata)'
    if r['contig'] == 'NODE_2578' and r['gene'] == 'amoC':
        r['note'] = 'comammox-type amoC (85.8% Ca. N. inopinata AmoC; 62.5% N. europaea); contig in bin.007 (DAS Tool)'
    if r['gene'] == 'nxrB' and r['contig'] in NXRB_NOTE and r['evidence'] == 'supported':
        r['note'] = '; '.join(x for x in (r['note'], NXRB_NOTE[r['contig']]) if x)
        if r['contig'] in NXRB_LIKE: r['gene'] = 'nxrB-like'

# not detected: unidades da Figura 1 sem acerto nem KOfam
UNITS = ['bin.004', 'bin.007', 'bin.009', 'bin.011', 'bin.016', 'bin.008', 'bin.017']
GENES_ND = ['amoA', 'amoB', 'amoC', 'nxrA/narG', 'hzsA', 'hzsB', 'hzsC', 'hdh', 'napA', 'nirS', 'nirK', 'nrfA', 'nirB', 'nirD', 'norB', 'norC', 'nosZ', 'mcrA', 'soxB', 'dsrA', 'dsrB']
present = {(r['bin'].rstrip('†‡'), 'nxrA/narG' if r['gene'] in ('narG', 'nxrA') else r['gene']) for r in rows}
for b in UNITS:
    for g in GENES_ND:
        if (b, g) in present: continue
        rows.append(dict(gene=g, bin=b, bin_taxon=genus(b), contig='', contig_length_bp='', depth_x='', contig_GC_pct='', hit_coordinates='', reference='',
                         identity_pct='', query_coverage_pct='', reference_coverage_pct='', evalue='', bit_score='', evidence='not detected', binner='',
                         DAS_Tool_retained='', DRAM_gene='', DRAM_KO='', note='no targeted hit >=40% identity and no KOfam call in this bin'))

ORDER = ['amoA', 'amoB', 'amoC', 'nxrA', 'narG', 'narG/nxrA (unassigned)', 'nxrA/narG', 'nxrB', 'nxrB-like', 'hzsA', 'hzsB', 'hzsC', 'hdh', 'napA', 'nirS', 'nirK', 'nrfA', 'nirB', 'nirD',
         'norB', 'norC', 'nosZ', 'mcrA', 'soxA', 'soxB', 'soxC', 'soxD', 'soxX', 'soxY', 'soxZ', 'dsrA', 'dsrB']
EVO = {'supported': 0, 'partial': 1, 'KOfam only': 2, 'not detected': 3}
df = pd.DataFrame(rows)
df['_g'] = df.gene.map({g: i for i, g in enumerate(ORDER)}); df['_e'] = df.evidence.map(EVO)
df['_b'] = df['bin'].map(lambda x: 'zzz' if x == 'unbinned' else x)
df = df.sort_values(['_g', '_e', '_b', 'identity_pct', 'contig'], ascending=[True, True, True, False, True], kind='mergesort').drop(columns=['_g', '_e', '_b'])
COLS = ['gene', 'bin', 'bin_taxon', 'contig', 'contig_length_bp', 'depth_x', 'contig_GC_pct', 'hit_coordinates', 'reference', 'identity_pct', 'tbn_identity_pct', 'orf_identity_pct', 'query_coverage_pct',
        'reference_coverage_pct', 'evalue', 'bit_score', 'evidence', 'binner', 'DAS_Tool_retained', 'DRAM_gene', 'DRAM_KO',
        'GC_z_vs_bin009', 'TNF_dist_to_bin009', 'TNF_null_P95_bin009', 'TNF_percentile_vs_bin009', 'note']
df = df.reindex(columns=COLS)
HDR = {'gene': 'Gene', 'bin': 'Bin', 'bin_taxon': 'Bin taxon (GTDB-Tk)', 'contig': 'Contig', 'contig_length_bp': 'Contig length (bp)', 'depth_x': 'Depth (x, samtools)',
       'contig_GC_pct': 'Contig GC (%)', 'hit_coordinates': 'Hit coordinates on contig', 'reference': 'Reference (accession)', 'identity_pct': 'Identity (%)', 'tbn_identity_pct': 'Identity tBLASTn vs assembly (%)', 'orf_identity_pct': 'Identity BLASTp on ORF (%)',
       'query_coverage_pct': 'Query coverage (%)', 'reference_coverage_pct': 'Reference coverage (%)', 'evalue': 'E-value', 'bit_score': 'Bit score',
       'evidence': 'Evidence', 'binner': 'Binner(s) that proposed the contig', 'DAS_Tool_retained': 'Retained by DAS Tool', 'DRAM_gene': 'DRAM gene at locus',
       'DRAM_KO': 'DRAM KO', 'GC_z_vs_bin009': f'GC z-score vs bin.009 ({GCM:.2f} ± {GCSD:.2f}%, unweighted, n=96)', 'TNF_dist_to_bin009': 'TNF distance to bin.009',
       'TNF_null_P95_bin009': 'TNF null P95 (bin.009 windows)', 'TNF_percentile_vs_bin009': 'TNF percentile vs bin.009 null', 'note': 'Note'}
out = df.rename(columns=HDR)
out.to_csv(os.path.join(HERE, 'TableS2_targeted_hits.tsv'), sep='\t', index=False)
pd.DataFrame(excluded).to_csv(os.path.join(HERE, 'TableS2_excluded_paralog_hits.tsv'), sep='\t', index=False)

# ---------------- aba References ----------------
import json
ORGJ = json.load(open(os.path.join(HERE, 'work', 'ref_org.json')))
ORGJ = {k.split('.')[0]: v for k, v in ORGJ.items()}
for a in list(SEQLEN):   # comprimento oficial do NCBI tem precedência
    if a.split('.')[0] in ORGJ and ORGJ[a.split('.')[0]]['slen']: SEQLEN[a] = int(ORGJ[a.split('.')[0]]['slen'])
SRC0 = lambda a: 'Swiss-Prot' if re.match(r'^[OPQ]\d|^[A-NR-Z]\d', a) and '_' not in a else ('RefSeq' if a.startswith(('WP_', 'NP_', 'YP_')) else 'GenBank (INSDC)')
def SRC(a):
    return 'PDB' if a == '8SR2_C' else SRC0(a)
def org(a):
    if a == 'Q1PZD8': return 'Candidatus Kuenenia stuttgartiensis (kustd1700)'
    if a.split('.')[0] in ORGJ:
        o = ORGJ[a.split('.')[0]]['organism']
        return o + (' (RefSeq MULTISPECIES; P. pantotrophus Sox cluster)' if a.startswith('WP_0248458') else '')
    m = re.search(r'\[([^\]]+)\]\s*$', DESC.get(a, '')); return m.group(1) if m else ''
REFS = []
use = {r: g for g, rs in GENES.items() for r in rs}
extra = {'WP_234989093.1': 'amoC (partial reference used in the original search)', 'WP_011421454.1': 'nosZ clade II', 'AAZ46324.1': 'nosZ clade II',
         'AAZ46320.1': 'nosZ clade II', 'WP_129545366.1': 'nosZ clade II', 'WP_301017512.1': 'nosZ clade II', 'GAB1413680.1': 'nosZ clade II',
         'ABI60535.1': 'amoA', 'WP_090322876.1': 'amoA', 'WP_011380154.1': 'amoA', 'WP_011330982.1': 'amoA', 'WP_062483313.1': 'amoC',
         'THJ16210.1': 'amoA', 'THJ16209.1': 'amoB', 'THJ18698.1': 'amoC', 'TKS61417.1': 'amoA', 'TKS61416.1': 'amoB', 'TKS61418.1': 'amoC',
         'Q82W83.1': 'amo gamma subunit (cytochrome c1)', 'WP_011111571.1': 'hao',
         'GAN34140.1': 'hzsA', 'GAN34138.1': 'hzsB', 'GAN34139.1': 'hzsC', 'WP_052564202.1': 'hdh', 'WP_080325085.1': 'hzsA', 'WP_080325083.1': 'hzsB',
         'WP_080325084.1': 'hzsC', 'OQD44929.1': 'hdh (HAO-like)', 'KKO20072.1': 'hzsA', 'KKO20886.1': 'hzsB (annotated)', 'KKO20885.1': 'hzsC (annotated)',
         'KKO18440.1': 'hdh (HAO-like)', 'P24474.1': 'nirS', 'WP_203415591.1': 'nirS', '8SR2_C': 'pmoC', 'WP_256036300.1': 'pmoC'}
for a, g in list(use.items()) + list(extra.items()):
    REFS.append(dict(accession=a, gene=g, organism=org(a), length_aa=SEQLEN.get(a, ''), source=SRC(a), description=re.sub(r'\s*\[[^\]]+\]\s*$', '', DESC.get(a, ''))[:90]))
# AmoA/PmoA da árvore de AmoA (Fig. S3; results/review_checks/comammox/amoA_tree/panel.tsv)
for l in open(os.path.join(RC, 'comammox', 'amoA_tree', 'panel.tsv')):
    a, grp, nm = l.rstrip('\n').split('\t')
    REFS.append(dict(accession=a, gene=f"{'pmoA' if grp == 'pmoA' else 'amoA'} (AmoA tree, Fig. S3; {grp.replace('_', ' ')})", organism=org(a) or nm,
                     length_aa=SEQLEN.get(a, '') or (ORGJ.get(a.split('.')[0], {}).get('slen') or ''), source=SRC(a), description=''))
# DsrA da árvore dsrA_typing
for l in open(os.path.join(RC, '11_sulfur', 'dsrA_tree_input.faa')):
    if l.startswith('>') and '|' in l:
        lab, a = l[1:].strip().split('|'); typ = {'OX': 'oxidative-type DsrA', 'RED': 'reductive-type DsrA (bacterial)', 'ARCH': 'reductive-type DsrA (archaeal; root)', 'TEST': 'DsrA (Ca. D. denitrificans)'}[lab.split('_')[0]]
        REFS.append(dict(accession=a, gene=f'dsrA ({typ}; dsrA_typing tree)', organism=org(a),
                         length_aa='', source=SRC(a), description=''))
seqlen_tree = {}
k = None
for l in open(os.path.join(RC, '11_sulfur', 'dsrA_tree_input.faa')):
    if l[0] == '>': k = l[1:].strip().split('|')[-1]; seqlen_tree[k] = 0
    elif k: seqlen_tree[k] += len(l.strip())
REF = pd.DataFrame(REFS)
REF['length_aa'] = [x if x != '' else seqlen_tree.get(a, '') for a, x in zip(REF.accession, REF.length_aa)]
REF.loc[REF.accession == 'Q1PZD8', 'organism'] = 'Candidatus Kuenenia stuttgartiensis (kustd1700)'
REF.loc[REF.accession == 'Q1PZD8', 'source'] = 'UniProtKB/TrEMBL'
REF = REF.drop_duplicates(['accession', 'gene'])
REF.columns = ['Accession', 'Gene / use', 'Organism', 'Length (aa)', 'Source', 'Description']
REF.to_csv(os.path.join(HERE, 'TableS2_references.tsv'), sep='\t', index=False)

# ---------------- XLSX ----------------
wb = Workbook(); ws = wb.active; ws.title = 'Table S2'
def sheet(ws, frame, widths):
    ws.append(list(frame.columns))
    for c in ws[1]: c.font = Font(bold=True); c.alignment = Alignment(wrap_text=True, vertical='top')
    for r in frame.itertuples(index=False): ws.append([None if (isinstance(x, float) and np.isnan(x)) or x == '' else x for x in r])
    for i, w in enumerate(widths): ws.column_dimensions[chr(65 + i) if i < 26 else 'A' + chr(65 + i - 26)].width = w
    ws.freeze_panes = 'B2'
sheet(ws, out, [10, 9, 16, 11, 10, 9, 8, 18, 15, 9, 9, 9, 9, 10, 9, 8, 12, 34, 9, 15, 9, 10, 10, 10, 10, 60])
for row in ws.iter_rows(min_row=2):
    for c in row:
        if isinstance(c.value, float): c.number_format = '0.0'
    for idx in (21, 22, 23):   # z e TNF com mais casas
        if isinstance(row[idx].value, float): row[idx].number_format = '0.00' if idx == 21 else '0.0000'
sheet(wb.create_sheet('References'), REF, [16, 34, 40, 10, 16, 60])
if excluded: sheet(wb.create_sheet('Excluded paralog hits'), pd.DataFrame(excluded), [8, 11, 9, 15, 9, 9, 9, 30, 90])
wb.save(os.path.join(HERE, 'TableS2_targeted_hits.xlsx'))
print(len(out), 'rows'); print(out['Evidence'].value_counts().to_string()); print(len(excluded), 'excluded paralog hits'); print(len(REF), 'references')
