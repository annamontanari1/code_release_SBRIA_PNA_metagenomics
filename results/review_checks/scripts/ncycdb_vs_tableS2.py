#!/usr/bin/env python3
"""NCycDB (95%) como terceira camada de comparacao com a busca dirigida (Table S2).

Etapas (cada uma pula se a saida ja existe; apague o arquivo para refazer):
  1. baixa o NCycDB (git clone de github.com/qichao1984/NCyc) e indexa NCyc_95 com DIAMOND
  2. prodigal -p meta na montagem completa (50.441 contigs; cobre bins e contigs sem bin)
  3. diamond blastp ORFs x NCyc_95 (E <= 1e-10, melhor hit por bitscore, inclui homologos nao-N)
  4. chamada NCycDB: melhor hit mapeado em id2gene.map, identidade >= 40%, cobertura da referencia >= 70%
     (mesmo criterio da busca dirigida); melhor hit fora do mapa = homologo (ex.: pmoA, hao-like, narG-like)
  5. concordancia com a Table S2: por linha (gene, contig, coordenadas), por bin e familias nao buscadas

Saidas: results/review_checks/14_ncycdb/  (SUMMARY.txt, TableS2_vs_NCycDB.tsv, compare_by_bin.tsv,
        ncyc_calls.tsv, ncyc_gene_vocabulary.tsv, ncyc_unsearched_families.tsv, excluded_paralogs_vs_NCycDB.tsv)
Ambiente: PATH com ~/miniconda3/envs/gunc/bin (diamond, prodigal, pandas); git do sistema.
"""
import os, sys, re, gzip, shutil, subprocess, time, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from base import *   # ROOT, RC, FINALBIN, node
import pandas as pd

THREADS = int(os.environ.get('THREADS', '8'))
DBDIR = os.path.expanduser('~/databases/NCyc')
OUT = os.environ.get('NCYC_OUT', RC + '/14_ncycdb').rstrip('/') + '/'; os.makedirs(OUT, exist_ok=True)
ASM = ROOT + '/results/targeted_search/bin_assemblies/NGS1140_6S_assembly.fa'
T2 = RC + '/supp_tables/TableS2_targeted_hits.tsv'
T2X = RC + '/supp_tables/TableS2_excluded_paralog_hits.tsv'
MIN_ID, MIN_SCOV, MAX_E = 40.0, 0.70, 1e-10
MIN_OVERLAP_NT = 60

def log(m): print(time.strftime('[%H:%M:%S] ') + m, flush=True)
def run(cmd, **kw):
    log('$ ' + ' '.join(map(str, cmd))); r = subprocess.run(cmd, **kw)
    if r.returncode: sys.exit(f'FALHOU ({r.returncode}): {cmd[0]}')

# ---------- 1. banco
def step_db():
    dmnd = DBDIR + '/NCyc_95.dmnd'
    if os.path.exists(dmnd) and os.path.exists(DBDIR + '/id2gene.map'): return dmnd
    if not os.path.isdir(DBDIR + '/.git'):
        os.makedirs(os.path.dirname(DBDIR), exist_ok=True)
        run(['git', 'clone', '--depth', '1', 'https://github.com/qichao1984/NCyc.git', DBDIR])
    found = {}
    for dp, _, fs in os.walk(DBDIR):
        for f in fs:
            if f in ('NCyc_95.faa', 'NCyc_95.faa.gz', 'id2gene.map'): found[f] = os.path.join(dp, f)
    log('arquivos do repositorio: ' + str(found))
    faa = found.get('NCyc_95.faa') or found.get('NCyc_95.faa.gz')
    if not faa or 'id2gene.map' not in found: sys.exit('NCyc_95.faa(.gz) ou id2gene.map nao encontrados no clone; veja ' + DBDIR)
    if os.path.getsize(faa) < 1_000_000: sys.exit(f'{faa} tem so {os.path.getsize(faa)} bytes (ponteiro LFS?). Baixe manualmente do repositorio.')
    plain = DBDIR + '/NCyc_95.faa'
    if faa.endswith('.gz'):
        with gzip.open(faa, 'rb') as a, open(plain, 'wb') as b: shutil.copyfileobj(a, b)
    elif faa != plain: shutil.copy(faa, plain)
    shutil.copy(found['id2gene.map'], DBDIR + '/id2gene.map')
    run(['diamond', 'makedb', '--in', plain, '-d', DBDIR + '/NCyc_95', '--threads', str(THREADS), '--quiet'])
    return dmnd

# ---------- 2/3. ORFs e busca
def _records(path):
    h, buf = None, []
    for l in open(path):
        if l[0] == '>':
            if h is not None: yield h, ''.join(buf)
            h, buf = l, []
        else: buf.append(l)
    if h is not None: yield h, ''.join(buf)

def _prodigal(fa, faa):
    r = subprocess.run(['prodigal', '-i', fa, '-a', faa, '-o', os.devnull, '-p', 'meta', '-q'], capture_output=True)
    return r.returncode

def step_orfs():
    """Prodigal em blocos de 1.000 contigs (em paralelo). O Prodigal x86 sob Rosetta falha (-11) com a montagem inteira.
    Blocos que falham sao refeitos contig a contig; contigs que ainda falham sao listados em orfs_failed_contigs.txt."""
    from concurrent.futures import ThreadPoolExecutor
    final = OUT + 'orfs_all.faa'
    if os.path.exists(final): return final
    cd = OUT + 'prodigal_chunks/'; os.makedirs(cd, exist_ok=True)
    chunks, cur, n = [], [], 0
    for rec in _records(ASM):
        cur.append(rec)
        if len(cur) == 1000: chunks.append(cur); cur = []
    if cur: chunks.append(cur)
    log(f'{sum(len(c) for c in chunks)} contigs em {len(chunks)} blocos')
    failed = []
    def work(i):
        fa, faa = f'{cd}c{i:03d}.fa', f'{cd}c{i:03d}.faa'
        if os.path.exists(faa + '.ok'): return
        with open(fa, 'w') as f:
            for h, sq in chunks[i]: f.write(h + sq)
        if _prodigal(fa, faa) == 0: open(faa + '.ok', 'w').close(); return
        log(f'bloco {i} falhou; refazendo contig a contig')
        with open(faa, 'w') as out:
            for k, (h, sq) in enumerate(chunks[i]):
                one, oneo = f'{cd}one_{i}.fa', f'{cd}one_{i}.faa'
                open(one, 'w').write(h + sq)
                if _prodigal(one, oneo) == 0: out.write(open(oneo).read())
                else: failed.append(h[1:].split()[0])
        open(faa + '.ok', 'w').close()
    with ThreadPoolExecutor(max_workers=max(1, THREADS)) as ex: list(ex.map(work, range(len(chunks))))
    with open(final + '.tmp', 'w') as out:
        for i in range(len(chunks)): out.write(open(f'{cd}c{i:03d}.faa').read())
    os.replace(final + '.tmp', final)
    open(OUT + 'orfs_failed_contigs.txt', 'w').write('\n'.join(failed) + ('\n' if failed else ''))
    log(f'ORFs gravadas em {final}; contigs sem predicao (Prodigal falhou): {len(failed)}')
    return final
def step_search(faa, dmnd):
    tsv = OUT + 'orfs_vs_NCyc95.tsv'
    if not os.path.exists(tsv):
        run(['diamond', 'blastp', '-q', faa, '-d', dmnd, '-o', tsv, '--more-sensitive', '-k', '10', '-e', str(MAX_E),
             '-p', str(THREADS), '--quiet', '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen',
             'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
    return tsv

# ---------- 4. chamadas
def load_orfs(faa):
    rows = []
    for l in open(faa):
        if l[0] != '>': continue
        p = l[1:].rstrip().split(' # ')
        rows.append((p[0].split()[0], node(p[0]), int(p[1]), int(p[2])))
    d = pd.DataFrame(rows, columns=['orf', 'contig', 'start', 'end'])
    d['bin'] = d.contig.map(lambda c: FINALBIN.get(c, 'unbinned')); return d

def load_map():
    m = {}
    for l in open(DBDIR + '/id2gene.map'):
        p = l.rstrip('\n').split('\t') if '\t' in l else l.split()
        if len(p) >= 2: m[p[0].strip()] = p[1].strip()
    return m
def canon(g): 
    x = re.match(r'[A-Za-z0-9]+', str(g)); return x.group(0).lower() if x else str(g).lower()

def calls(tsv, orfs, idmap):
    cols = ['orf', 'sseqid', 'pident', 'length', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    b = pd.read_csv(tsv, sep='\t', names=cols)
    b = b.sort_values('bitscore', ascending=False).drop_duplicates('orf')
    ids = set(b.sseqid)
    log(f'IDs de subject no melhor hit: {len(ids)}; presentes em id2gene.map: {len(ids & set(idmap))}; '
        f'exemplos subject={sorted(ids)[:3]} map={list(idmap)[:3]}')
    if idmap and not (ids & set(idmap)): sys.exit('Nenhum ID do banco casa com id2gene.map: o formato dos IDs difere. Veja os exemplos acima.')
    b['scov'] = ((b.send - b.sstart).abs() + 1) / b.slen
    b['gene'] = b.sseqid.map(idmap)
    b['canon'] = b.gene.map(lambda g: canon(g) if isinstance(g, str) else None)
    def cls(r):
        if not isinstance(r.gene, str): return 'homolog_not_in_map'
        return 'ncyc_call' if (r.pident >= MIN_ID and r.scov >= MIN_SCOV) else 'ncyc_below_threshold'
    b['cls'] = b.apply(cls, axis=1)
    b = b.merge(orfs, on='orf', how='left'); return b

# ---------- 5. comparacao
def tokens(g):
    t = {x.lower() for x in re.findall(r'[A-Za-z][A-Za-z0-9]*', str(g))}; return t - {'unassigned'}

# explicacao das discordancias (21/10/2026; 14_ncycdb/followup_21out). Sem categoria -> '' e lista a parte (discordance_unexplained.tsv)
OUT_OF_SCOPE = ('sox', 'dsr', 'mcr')
EXPLAINED = ('discordant_other_gene', 'homolog_not_in_map', 'other_gene_below_threshold')
HOMOLOG_SAME = 'homolog of the same family (not in NCycDB gene map)'
def same_family(tk, n):
    return n in tk or bool(tk & {'narg', 'nxra'} and n in ('narz', 'narg', 'nxra')) or ('nxrb' in tk and n in ('narh', 'nary', 'nxrb'))
def explain(gene, ng, result, hom_desc=None):
    # hom_desc: gene da descricao do cabecalho do NCyc_95.faa ("[description=narG homolog ...]") quando o melhor hit esta fora do mapa
    if result not in EXPLAINED: return ''
    tk = tokens(gene); n = canon(ng) if isinstance(ng, str) and not ng.startswith('(') else None
    if any(t.startswith(OUT_OF_SCOPE) for t in tk): return 'outside NCycDB scope'
    if result == 'homolog_not_in_map': return HOMOLOG_SAME if hom_desc and same_family(tk, canon(hom_desc)) else ''
    if tk & {'amoa', 'amob', 'amoc'} and n in ('pmoa', 'pmob', 'pmoc'): return 'pmo'
    if tk & {'narg', 'nxra'} and n in ('narz', 'narg', 'nxra'): return 'nar/nxr family'
    if 'nxrb' in tk and n in ('narh', 'nary', 'nxrb'): return 'nar/nxr family'   # inclui NxrB de anammox dentro de 'narH' do NCycDB
    if tk & {'hzsb', 'hzsc'} and n in ('hzsb', 'hzsc') and n not in tk: return 'hzsB/hzsC nomenclature'
    return ''

def compare(b, orfs):
    t2 = pd.read_csv(T2, sep='\t'); t2 = t2[t2['Gene'].notna()]
    by_contig = {c: g for c, g in b.groupby('contig')}
    res = []
    for _, r in t2.iterrows():
        tk = tokens(r['Gene']); m = re.match(r'(\d+)-(\d+)', str(r['Hit coordinates on contig']))
        if pd.isna(r['Contig']) or str(r['Contig']).strip() == '':   # linha sem contig (ex.: not detected): confere no nivel do bin
            bn = str(r['Bin']).replace('‡', '')
            n = int(((b.cls == 'ncyc_call') & b.canon.isin(tk) & (b['bin'] == bn)).sum())
            res.append(dict(Gene=r['Gene'], Bin=bn, Contig='', Evidence=r.get('Evidence', ''), TableS2_identity='',
                            NCycDB_best=f'{n} ORF(s) no bin' if n else '', NCycDB_identity='', NCycDB_ref_cov='',
                            result='bin_level_NCycDB_calls_gene' if n else 'bin_level_consistent_no_call'))
            continue
        s, e = (int(m.group(1)), int(m.group(2))) if m else (None, None)
        g = by_contig.get(r['Contig'])
        sub = None
        if g is not None and s is not None:
            ov = g[(g[['end']].values.ravel().clip(max=e) - g[['start']].values.ravel().clip(min=s)) >= MIN_OVERLAP_NT]
            sub = ov.sort_values('bitscore', ascending=False)
        sid = ''
        if sub is None or sub.empty: out, ng, pid, sc = 'no_NCycDB_hit', '', '', ''
        else:
            top = sub.iloc[0]; ng = top.gene if isinstance(top.gene, str) else '(homolog fora do mapa)'
            pid, sc = round(top.pident, 1), round(top.scov * 100, 1)
            same = isinstance(top.gene, str) and top.canon in tk
            out = ('concordant' if same and top.cls == 'ncyc_call' else
                   'same_gene_below_threshold' if same else
                   'homolog_not_in_map' if top.cls == 'homolog_not_in_map' else
                   'discordant_other_gene' if top.cls == 'ncyc_call' else 'other_gene_below_threshold')
            sid = top.sseqid
        res.append(dict(Gene=r['Gene'], Bin=str(r['Bin']).replace('‡', ''), Contig=r['Contig'], Evidence=r.get('Evidence', ''),
                        TableS2_identity=r.get('Identity BLASTp on ORF (%)', ''), NCycDB_best=ng, NCycDB_identity=pid,
                        NCycDB_ref_cov=sc, result=out, _sseqid=sid))
    R = pd.DataFrame(res)
    # cabecalhos do NCyc_95.faa para os subjects das linhas que podem ser explicadas (descricao dos homologos fora do mapa)
    want, hdr = set(R[R.result.isin(EXPLAINED)]._sseqid), {}
    for l in open(DBDIR + '/NCyc_95.faa'):
        if l[0] == '>':
            p = l[1:].rstrip().split(' ', 1)
            if p[0] in want: hdr[p[0]] = p[1] if len(p) > 1 else ''
    hd = lambda x: (re.search(r'description=(\S+) homolog', hdr.get(x, '')) or [None, None])[1]
    R['discordance_explanation'] = [explain(g, ng, o, hd(x)) for g, ng, o, x in zip(R.Gene, R.NCycDB_best, R.result, R._sseqid)]
    R.drop(columns='_sseqid').to_csv(OUT + 'TableS2_vs_NCycDB.tsv', sep='\t', index=False)
    # linhas sem categoria: listadas a parte, com o subject do NCycDB e a descricao do cabecalho do NCyc_95.faa
    U = R[R.result.isin(EXPLAINED) & (R.discordance_explanation == '')].copy()
    U['NCycDB_subject'] = U._sseqid; U['NCycDB_subject_header'] = U._sseqid.map(lambda x: hdr.get(x, ''))
    U.drop(columns=['_sseqid', 'discordance_explanation']).to_csv(OUT + 'discordance_unexplained.tsv', sep='\t', index=False)
    R = R.drop(columns='_sseqid')
    # excluidos como parálogos: o NCycDB os chamaria como o gene-alvo?
    x = pd.read_csv(T2X, sep='\t'); xr = []
    for _, r in x.iterrows():
        tk = tokens(r['gene']); g = by_contig.get(r['contig'])
        hit = g[(g.cls == 'ncyc_call') & (g.canon.isin(tk))] if g is not None else None
        xr.append(dict(gene=r['gene'], contig=r['contig'], bin=r['bin'], reason=r['reason'],
                       NCycDB_calls_target_gene='yes' if hit is not None and len(hit) else 'no',
                       NCycDB_genes_on_contig=';'.join(sorted(set(g[g.cls == 'ncyc_call'].gene))) if g is not None else ''))
    X = pd.DataFrame(xr); X.to_csv(OUT + 'excluded_paralogs_vs_NCycDB.tsv', sep='\t', index=False)
    # por bin
    c = b[b.cls == 'ncyc_call']
    nb = c.groupby(['canon', 'bin']).size().rename('n_ORFs_NCycDB').reset_index()
    R['canon'] = R.Gene.map(lambda g: sorted(tokens(g))[0] if tokens(g) else '')
    tb = R.assign(canon=R.Gene.map(lambda g: '|'.join(sorted(tokens(g))))).groupby(['canon', 'Bin']).size().rename('n_TableS2_rows').reset_index()
    tb = tb.rename(columns={'Bin': 'bin'})
    tb['bin'] = tb['bin'].replace({'nan': 'unbinned', '': 'unbinned'})
    nb['canon'] = nb.canon.astype(str)
    cmp_ = nb.merge(tb, on=['canon', 'bin'], how='outer').fillna(0); cmp_.to_csv(OUT + 'compare_by_bin.tsv', sep='\t', index=False)
    allt2 = set().union(*[tokens(g) for g in t2['Gene']])
    un = c[~c.canon.isin(allt2)].groupby(['gene', 'bin']).size().rename('n_ORFs').reset_index()
    un.to_csv(OUT + 'ncyc_unsearched_families.tsv', sep='\t', index=False)
    return R, X, un

def main():
    dmnd = step_db(); faa = step_orfs(); tsv = step_search(faa, dmnd)
    orfs = load_orfs(faa); idmap = load_map(); b = calls(tsv, orfs, idmap)
    b[b.cls == 'ncyc_call'][['orf', 'contig', 'bin', 'gene', 'sseqid', 'pident', 'scov', 'evalue', 'bitscore']].to_csv(OUT + 'ncyc_calls.tsv', sep='\t', index=False)
    pd.Series(collections.Counter(idmap.values())).sort_index().rename('n_ref_seqs').to_csv(OUT + 'ncyc_gene_vocabulary.tsv', sep='\t')
    R, X, un = compare(b, orfs)
    with open(OUT + 'SUMMARY.txt', 'w') as f:
        f.write(f'ORFs: {len(orfs)}; ORFs com hit no NCyc_95: {len(b)}; chamadas NCycDB (id>={MIN_ID}%, cob. ref>={MIN_SCOV*100:.0f}%): {(b.cls=="ncyc_call").sum()}; '
                f'melhor hit fora do mapa (homologos): {(b.cls=="homolog_not_in_map").sum()}\n\n')
        f.write('Table S2 x NCycDB (linhas da Table S2 por resultado):\n' + R.result.value_counts().to_string() + '\n\n')
        f.write('Por evidencia da Table S2:\n' + pd.crosstab(R.Evidence, R.result).to_string() + '\n\n')
        f.write('Por gene (Table S2):\n' + pd.crosstab(R.Gene, R.result).to_string() + '\n\n')
        d = R[R.result.isin(EXPLAINED)]
        f.write('Explicacao (discordant_other_gene, homolog_not_in_map, other_gene_below_threshold; vazio = sem categoria, ver discordance_unexplained.tsv):\n'
                + pd.crosstab(d.discordance_explanation.replace('', '(sem categoria)'), d.result).to_string() + '\n\n')
        s = R[R.Evidence == 'supported']; scope = s.Gene.map(lambda g: any(t.startswith(OUT_OF_SCOPE) for t in tokens(g)))
        si = s[~scope]
        f.write(f'Linhas supported: {len(s)}; fora do escopo do NCycDB (sox, dsr, mcr): {int(scope.sum())} '
                f'({s[scope].result.value_counts().to_dict()}); restantes: {len(si)}\n')
        dd = si[si.result == 'discordant_other_gene']
        f.write(f"  concordant: {(si.result=='concordant').sum()}\n  same_gene_below_threshold: {(si.result=='same_gene_below_threshold').sum()}\n"
                f"  discordant explicada: {(dd.discordance_explanation!='').sum()} {dd[dd.discordance_explanation!=''].discordance_explanation.value_counts().to_dict()}\n"
                f"  discordant nao explicada: {(dd.discordance_explanation=='').sum()}\n"
                f"  homolog_not_in_map: {(si.result=='homolog_not_in_map').sum()} (explicadas: {((si.result=='homolog_not_in_map') & (si.discordance_explanation!='')).sum()})\n"
                f"  other_gene_below_threshold: {(si.result=='other_gene_below_threshold').sum()} (explicadas: {((si.result=='other_gene_below_threshold') & (si.discordance_explanation!='')).sum()})\n"
                f"  no_NCycDB_hit: {(si.result=='no_NCycDB_hit').sum()}\n")
        f.write('Supported no escopo (' + str(len(si)) + ') por resultado x categoria:\n'
                + pd.crosstab(si.result, si.discordance_explanation.replace('', '-')).to_string() + '\n'
                + f"  sem explicacao (resultado discordant/homolog/other_below sem categoria): {(si.result.isin(EXPLAINED) & (si.discordance_explanation=='')).sum()}\n\n")
        f.write(f"Paralogos excluidos que o NCycDB chamaria como o gene-alvo: {(X.NCycDB_calls_target_gene=='yes').sum()} de {len(X)}\n\n")
        f.write('Familias do NCycDB nao buscadas na Table S2 (ORFs chamadas por gene):\n' + (un.groupby('gene').n_ORFs.sum().sort_values(ascending=False).to_string() if len(un) else '(nenhuma)') + '\n')
    print(open(OUT + 'SUMMARY.txt').read()); log('Saidas em ' + OUT)

if __name__ == '__main__': main()
