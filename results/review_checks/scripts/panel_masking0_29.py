#!/usr/bin/env python3
"""Seção 29 do README (24/10/2026): DIAMOND do painel anammox de 22 sequências (seção 12) sem mascaramento e com rótulos corrigidos.
Não altera 12_anammox_contigs/ nem a seção 12. Saídas em 29_panel_masking0/.
Mesmas entradas e parâmetros de step_ref_identity() em classify_anammox_contigs.py (DIAMOND 2.1.24 do ambiente gunc; ORFs orfs10.faa;
painel anammox_refs.faa; --more-sensitive -k 200 -e 1e-10; pident >= 40; subunidade = rótulo da referência de maior bit score;
melhor hit por gênero = maior pident dentro da subunidade). Variam apenas --masking e --comp-based-stats (padrões da 2.1.24: tantan e 1).
Rótulos corrigidos (PDB 5C2V / Swiss-Prot): KKO20885.1 = hzsB (beta), KKO20886.1 = hzsC (gamma); o GenBank traz o inverso.
Uso (de ~/sbria_anammox_shotgun): python results/review_checks/scripts/panel_masking0_29.py
"""
import os, re, subprocess
import pandas as pd

RC = 'results/review_checks'
S12 = f'{RC}/12_anammox_contigs'
OUT = f'{RC}/29_panel_masking0'
DIAMOND = os.path.expanduser('~/miniconda3/envs/gunc/bin/diamond')
ORFS, REFS = f'{S12}/orfs10.faa', f'{S12}/anammox_refs.faa'
FIX = {'KKO20885.1': 'hzsB', 'KKO20886.1': 'hzsC'}
COLS = ['orf', 'acc', 'pident', 'alnlen', 'qlen', 'slen', 'qs', 'qe', 'ss', 'se', 'evalue', 'bitscore']
VARIANTS = {   # nome: (masking, comp-based-stats); m1c1 = padrão (= seção 12)
    'm1c1': ('tantan', '1'), 'm0c1': ('0', '1'), 'm1c0': ('tantan', '0'), 'm0c0': ('0', '0')}
MAIN = 'm0c0'  # --masking 0 e --comp-based-stats 0 (a 2.1.24 aplica ajuste composicional por padrão)
PUBLISHED = {'hzsA': (88.2, 97.5), 'hzsB': (94.3, 96.3), 'hzsC': (78.7, 98.2), 'hdh': (94.8, 97.7)}
PUB_MARGIN = (5.9, 17.5)
os.makedirs(OUT, exist_ok=True)


def run_diamond(name, masking, cbs):
    tsv = f'{OUT}/diamond_{name}.tsv'
    subprocess.run([DIAMOND, 'blastp', '-q', ORFS, '-d', f'{OUT}/anammox_refs', '-o', tsv, '--more-sensitive', '-k', '200',
                    '-e', '1e-10', '-p', '4', '--quiet', '--masking', masking, '--comp-based-stats', cbs,
                    '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send',
                    'evalue', 'bitscore'], check=True)
    return pd.read_csv(tsv, sep='\t', header=None, names=COLS)


def by_genus(h, tab, labels):
    """Mesma lógica de step_ref_identity(); labels = coluna de subunidade usada ('gene_genbank' ou 'gene_fixed')."""
    h = h.merge(tab[['accession', labels, 'genus']].rename(columns={labels: 'gene'}), left_on='acc', right_on='accession')
    h = h[h.pident >= 40]
    rows = []
    for orf, g in h.groupby('orf'):
        top = g.sort_values('bitscore', ascending=False).iloc[0]
        gg = g[g.gene == top.gene]
        per = gg.sort_values('pident', ascending=False).groupby('genus').first().sort_values('pident', ascending=False)
        br = per.loc['Brocadia'] if 'Brocadia' in per.index else None
        ku = per.loc['Kuenenia'] if 'Kuenenia' in per.index else None
        rows.append(dict(orf=orf, gene=top.gene, top_bitscore_ref=top.acc, orf_len_aa=int(g.qlen.iloc[0]),
                         best_genus=per.index[0],
                         brocadia_ref=br.acc if br is not None else '', brocadia_pident=round(br.pident, 1) if br is not None else None,
                         kuenenia_ref=ku.acc if ku is not None else '', kuenenia_pident=round(ku.pident, 1) if ku is not None else None,
                         margin_vs_kuenenia=round(br.pident - ku.pident, 1) if br is not None and ku is not None else None,
                         all_genera=';'.join('%s:%.1f' % (k, v.pident) for k, v in per.iterrows())))
    return pd.DataFrame(rows)


def ranges(df):
    d = df[df.orf_len_aa >= 100]
    out = {}
    for gene in ['hzsA', 'hzsB', 'hzsC', 'hdh']:
        s = d[d.gene == gene]
        out[gene] = (s.brocadia_pident.min(), s.brocadia_pident.max(), len(s))
    m = d.margin_vs_kuenenia.dropna()
    return out, (m.min(), m.max(), len(m), int(d.margin_vs_kuenenia.isna().sum()))


def main():
    subprocess.run([DIAMOND, 'makedb', '--in', REFS, '-d', f'{OUT}/anammox_refs', '--quiet'], check=True)
    tab = pd.read_csv(f'{S12}/anammox_reference_panel.tsv', sep='\t').rename(columns={'gene': 'gene_genbank'})
    tab['gene_fixed'] = [FIX.get(a, g) for a, g in zip(tab.accession, tab.gene_genbank)]
    tab.to_csv(f'{OUT}/panel_labels_corrected.tsv', sep='\t', index=False)
    hits = {n: run_diamond(n, *v) for n, v in VARIANTS.items()}

    # controle: a variante padrão reproduz a saída da seção 12
    old = pd.read_csv(f'{S12}/orfs_vs_anammox_refs.tsv', sep='\t', header=None).iloc[:, :12]
    old.columns = COLS
    key = ['orf', 'acc', 'pident', 'alnlen', 'qs', 'qe', 'ss', 'se']
    same = old[key].sort_values(key).reset_index(drop=True).equals(hits['m1c1'][key].sort_values(key).reset_index(drop=True))

    # subunidade pela sequência: cadeia de 5C2V de maior bit score (BLASTp existente, followup_21out) + ordem no operon
    c5 = pd.read_csv(f'{RC}/14_ncycdb/followup_21out/hzs_orfs_vs_hzsBC_refs.tsv', sep='\t', header=None)
    c5 = c5[c5[1].str.contains(r'Q1Q0T[34]')].sort_values(8, ascending=False).drop_duplicates(0)
    chain = {o: ('hzsB' if 'Q1Q0T4' in s else 'hzsC', p) for o, s, p in zip(c5[0], c5[1], c5[2])}
    gff = []
    for l in open(f'{S12}/prodigal.gff'):
        if l.startswith('#') or '\tCDS\t' not in l:
            continue
        f = l.rstrip('\n').split('\t')
        n = re.search(r'ID=\d+_(\d+)', f[8]).group(1)
        part = re.search(r'partial=(\d\d)', f[8]).group(1)
        gff.append(dict(orf=f'{f[0]}_{n}', contig=f[0], start=int(f[3]), end=int(f[4]), strand=f[6], prodigal_partial=part))
    gff = pd.DataFrame(gff)

    corr = {n: by_genus(h, tab, 'gene_fixed') for n, h in hits.items()}
    pub = by_genus(hits['m1c1'], tab, 'gene_genbank')
    ref_pub = pd.read_csv(f'{S12}/orf_identity_by_genus.tsv', sep='\t')
    reproduces_pub = (pub.set_index('orf').loc[ref_pub.orf, 'brocadia_pident'].values == ref_pub.best_pident.values).all()

    a = corr[MAIN].merge(gff, on='orf', how='left')
    a['label_sec12'] = a.orf.map(pub.set_index('orf').gene)
    a['chain_5C2V'] = a.orf.map(lambda o: chain.get(o, ('', None))[0])
    a['chain_5C2V_pident'] = a.orf.map(lambda o: chain.get(o, ('', None))[1])
    pos = {}
    for c, g in a[a.gene.str.startswith('hzs')].groupby('contig'):
        g = g.sort_values('start', ascending=(g.strand.iloc[0] == '+'))
        for i, o in enumerate(g.orf, 1):
            pos[o] = i
        pos[c] = '-'.join(g.gene.str[-1])
    a['operon_position'] = a.orf.map(pos)
    a['operon_order'] = a.contig.map(lambda c: pos.get(c, ''))
    a['subunit_by_sequence'] = [ch if ch else g for ch, g in zip(a.chain_5C2V, a.gene)]
    a['consistent'] = a.subunit_by_sequence == a.gene
    a['label_changed'] = a.gene != a.label_sec12
    a['partial'] = a.prodigal_partial.map(lambda p: 'no' if p == '00' else 'yes')
    a = a.sort_values(['gene', 'contig', 'orf'])
    a.to_csv(f'{OUT}/orf_assignment_masking0.tsv', sep='\t', index=False)
    for n, d in corr.items():
        d.to_csv(f'{OUT}/orf_identity_by_genus_{n}_labels_corrected.tsv', sep='\t', index=False)

    # 60 pares ORF-referência presentes no DIAMOND e no BLASTp do 06_tableS2
    bl = pd.read_csv(f'{RC}/06_tableS2/hzs_orfs_vs_anammox_panel_long.tsv', sep='\t', header=None,
                     names=['orf', 'acc', 'pident', 'alnlen', 'qs', 'qe', 'ss', 'se', 'evalue', 'bitscore', 'qlen', 'slen'])
    best = lambda h: h.sort_values('bitscore', ascending=False).drop_duplicates(['orf', 'acc']).set_index(['orf', 'acc']).pident
    p = pd.DataFrame({'sec12_published': best(old), 'default_labels_corrected': best(hits['m1c1']),
                      'masking0_cbs1': best(hits['m0c1']), 'masking_default_cbs0': best(hits['m1c0']),
                      'masking0_cbs0': best(hits['m0c0']), 'blastp_06_tableS2': best(bl)}).dropna(subset=['sec12_published', 'blastp_06_tableS2'])
    p = p.round(1).reset_index()
    p['gene_fixed'] = p.acc.map(tab.set_index('accession').gene_fixed)
    p['d_masking0_vs_default'] = (p.masking0_cbs0 - p.default_labels_corrected).round(1)
    p['d_masking0_vs_blastp'] = (p.masking0_cbs0 - p.blastp_06_tableS2).round(1)
    p['changes_gt_0.5'] = p.d_masking0_vs_default.abs() > 0.5
    p.to_csv(f'{OUT}/pairs60_identity_by_variant.tsv', sep='\t', index=False)

    L = [f'DIAMOND {subprocess.run([DIAMOND, "version"], capture_output=True, text=True).stdout.strip()}; variante principal {MAIN} '
         f'(--masking {VARIANTS[MAIN][0]} --comp-based-stats {VARIANTS[MAIN][1]}); padrão = --masking tantan --comp-based-stats 1',
         f'controle: variante padrão reproduz orfs_vs_anammox_refs.tsv da seção 12: {same}; '
         f'rótulos do GenBank + padrão reproduzem orf_identity_by_genus.tsv (best_pident): {bool(reproduces_pub)}',
         f'pares ORF-referência comparáveis: {len(p)}; mudam > 0,5 ponto (masking0_cbs0 vs padrão): {int(p["changes_gt_0.5"].sum())}; '
         f'máx |masking0_cbs0 - BLASTp| = {p.d_masking0_vs_blastp.abs().max()}', '']
    r_pub, m_pub = ranges(pub)
    L.append('faixas (ORFs >= 100 aa; melhor Brocadia) e margem sobre a melhor Kuenenia da mesma subunidade:')
    for n, label in [('pub', 'seção 12 (rótulos GenBank, padrão)'), ('m1c1', 'rótulos corrigidos, padrão'), ('m0c0', 'rótulos corrigidos, masking 0 / cbs 0'),
                     ('m0c1', 'rótulos corrigidos, masking 0 / cbs 1')]:
        r, m = ranges(pub if n == 'pub' else corr[n])
        L.append(f'  {label}: ' + '; '.join(f'{g} {v[0]:.1f}-{v[1]:.1f} (n={v[2]})' for g, v in r.items())
                 + f'; margem {m[0]:.1f}-{m[1]:.1f} (n={m[2]}; sem Kuenenia: {m[3]})')
    L.append('  publicado: ' + '; '.join(f'{g} {v[0]}-{v[1]}' for g, v in PUBLISHED.items()) + f'; margem {PUB_MARGIN[0]}-{PUB_MARGIN[1]}')
    open(f'{OUT}/SUMMARY.txt', 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
