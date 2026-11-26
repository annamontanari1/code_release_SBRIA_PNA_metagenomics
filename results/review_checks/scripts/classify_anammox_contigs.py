#!/usr/bin/env python3
"""
Classifica taxonomicamente os 10 contigs com hzsA/hzsB/hzsC/hdh (nao atribuidos ao bin.009).

Etapas
 1. extrai os 10 contigs do assembly
 2. preve ORFs (prodigal -p meta)
 3. melhor hit taxonomico de cada ORF no proGenomes 2.1 (banco do GUNC; DIAMOND --more-sensitive)
 4. monta painel de referencias hzsA/B/C e hdh de Brocadia, Kuenenia, Jettenia, Scalindua, Anammoxoglobus
    (referencias locais + NCBI E-utilities) e mede a identidade de cada ORF por genero
 5. (opcional, se mafft/trimal/iqtree estiverem no PATH) arvores de hzsA e hdh
Saidas em results/review_checks/12_anammox_contigs/
Uso:  export PATH=~/miniconda3/envs/gunc/bin:$PATH ; python classify_anammox_contigs.py
"""
import os, re, sys, glob, shutil, subprocess, time, urllib.parse, urllib.request
import pandas as pd

P = os.path.expanduser('~/sbria_anammox_shotgun')
ASM = f'{P}/results/targeted_search/bin_assemblies/NGS1140_6S_assembly.fa'
LOCAL_REFS = f'{P}/results/review_checks/refs/refs_clean.faa'
OUT = f'{P}/results/review_checks/12_anammox_contigs'
THREADS = int(os.environ.get('THREADS', '8'))
CONTIGS = ['NODE_5668', 'NODE_5312', 'NODE_7170', 'NODE_17695', 'NODE_29475',
           'NODE_24575', 'NODE_20914', 'NODE_23806', 'NODE_47442', 'NODE_24842']
GENERA = ['Brocadia', 'Kuenenia', 'Jettenia', 'Scalindua', 'Anammoxoglobus', 'Anammoximicrobium']
GENE_TERMS = {'hzsA': ['hydrazine synthase subunit alpha', 'hydrazine synthase alpha subunit'],
              'hzsB': ['hydrazine synthase subunit beta', 'hydrazine synthase beta subunit'],
              'hzsC': ['hydrazine synthase subunit gamma', 'hydrazine synthase gamma subunit'],
              'hdh': ['hydrazine dehydrogenase']}
KNOWN_GENUS = {'Q1PW30.2': 'Kuenenia'}  # HDH de K. stuttgartiensis (Table S2)
os.makedirs(OUT, exist_ok=True)


def log(*a):
    print('[%s]' % time.strftime('%H:%M:%S'), *a, flush=True)


def need(tool):
    p = shutil.which(tool)
    if not p:
        sys.exit(f'ERRO: "{tool}" nao encontrado no PATH. Rode antes: export PATH=~/miniconda3/envs/gunc/bin:$PATH')
    return p


def read_fasta(path):
    name, seq = None, []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(seq)
                name, seq = line[1:].strip(), []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, ''.join(seq)


def write_fasta(path, records):
    with open(path, 'w') as f:
        for h, s in records:
            f.write('>%s\n' % h)
            for i in range(0, len(s), 80):
                f.write(s[i:i + 80] + '\n')


def run(cmd, **kw):
    log('$', ' '.join(map(str, cmd)))
    subprocess.run(cmd, check=True, **kw)


# ---------------------------------------------------------------- 1. contigs
def step_contigs():
    out = f'{OUT}/contigs10.fna'
    want, recs = set(CONTIGS), []
    for h, s in read_fasta(ASM):
        node = h.split()[0].split('_length')[0]
        if node in want:
            recs.append((node, s))
    got = {r[0] for r in recs}
    if got != want:
        sys.exit('ERRO: contigs ausentes no assembly: %s' % sorted(want - got))
    write_fasta(out, recs)
    log('contigs extraidos:', ', '.join('%s (%d bp)' % (n, len(s)) for n, s in recs))
    return out


# ---------------------------------------------------------------- 2. ORFs
def step_orfs(fna):
    faa = f'{OUT}/orfs10.faa'
    run([need('prodigal'), '-i', fna, '-a', faa, '-p', 'meta', '-q', '-o', f'{OUT}/prodigal.gff', '-f', 'gff'])
    recs = [(h.split()[0], s.rstrip('*')) for h, s in read_fasta(faa)]
    write_fasta(faa, recs)
    log('ORFs previstas:', len(recs))
    return faa


# ---------------------------------------------------------------- 3. proGenomes
def find_one(patterns):
    for pat in patterns:
        hits = glob.glob(os.path.expanduser(pat), recursive=True)
        if hits:
            return hits[0]
    return None


def step_progenomes(orfs):
    db = find_one(['~/databases/**/gunc_db_progenomes2.1.dmnd', '~/databases/gunc_db_progenomes2.1.dmnd'])
    tax = find_one(['~/databases/**/genome2taxonomy_pg2.1ref.tsv',
                    '~/miniconda3/envs/gunc/lib/python*/site-packages/gunc/**/genome2taxonomy_pg2.1ref.tsv'])
    if not db or not tax:
        sys.exit('ERRO: nao achei o banco do GUNC (%s) ou o mapa de taxonomia (%s). Procure com:\n'
                 '  find ~ -name "genome2taxonomy_pg2.1ref.tsv" 2>/dev/null' % (db, tax))
    tsv = f'{OUT}/orfs_vs_progenomes2.1.tsv'
    if not os.path.exists(tsv) or os.path.getsize(tsv) == 0:
        run([need('diamond'), 'blastp', '-q', orfs, '-d', db, '-o', tsv, '--more-sensitive', '-k', '5', '-e', '1e-5',
             '-b', '1.0', '-c', '4', '-p', str(THREADS), '--tmpdir', OUT,
             '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen', 'evalue', 'bitscore'])
    hits = pd.read_csv(tsv, sep='\t', header=None,
                       names=['orf', 'sseqid', 'pident', 'alnlen', 'qlen', 'slen', 'evalue', 'bitscore'])
    m = pd.read_csv(tax, sep='\t', dtype=str)
    log('mapa de taxonomia:', tax, '| colunas:', list(m.columns))
    idcol = m.columns[0]
    cols = {c.lower(): c for c in m.columns}
    ranks = {r: cols[r] for r in ['phylum', 'class', 'order', 'family', 'genus', 'species'] if r in cols}
    if 'genus' not in ranks:
        lin = m.columns[-1]
        def pull(s, p):
            mm = re.search(r'%s__([^;|]*)' % p, str(s))
            return mm.group(1) if mm else ''
        for r, p in [('phylum', 'p'), ('family', 'f'), ('genus', 'g'), ('species', 's')]:
            m[r] = m[lin].map(lambda s, p=p: pull(s, p)); ranks[r] = r
        if (m['genus'] == '').all():
            print(m.head(3).to_string()); sys.exit('ERRO: formato do mapa de taxonomia nao reconhecido (cole as 3 linhas acima).')
    mm = m.set_index(idcol)
    def lookup(s):
        for key in (s, s.rsplit('_', 1)[0], s.split('_')[0]):
            if key in mm.index:
                r = mm.loc[key]
                r = r.iloc[0] if isinstance(r, pd.DataFrame) else r
                return {k: r[v] for k, v in ranks.items()}
        return {k: '' for k in ranks}
    hits = hits.sort_values(['orf', 'bitscore'], ascending=[True, False])
    rows = []
    for orf, g in hits.groupby('orf', sort=False):
        top = g.iloc[0]
        t = lookup(top.sseqid)
        gen = [lookup(x).get('genus', '') for x in g.sseqid]
        rows.append(dict(orf=orf, contig=orf.rsplit('_', 1)[0], best_genome=top.sseqid, pident=top.pident,
                         qcov=round(top.alnlen / top.qlen, 2), bitscore=top.bitscore,
                         phylum=t.get('phylum', ''), family=t.get('family', ''), genus=t.get('genus', ''),
                         species=t.get('species', ''), top5_genera=';'.join(x or '?' for x in gen)))
    df = pd.DataFrame(rows)
    df.to_csv(f'{OUT}/orf_besthit_progenomes.tsv', sep='\t', index=False)
    return df


# ---------------------------------------------------------------- 4. referencias anammox
def classify_header(h):
    low = h.lower()
    gene = None
    if 'hydrazine dehydrogenase' in low or re.search(r'\bhdh\b', low):
        gene = 'hdh'
    elif 'hydrazine synthase' in low or 'hzs' in low:
        if 'alpha' in low or re.search(r'subunit\s*a\b', low):
            gene = 'hzsA'
        elif 'beta' in low or re.search(r'subunit\s*b\b', low):
            gene = 'hzsB'
        elif 'gamma' in low or re.search(r'subunit\s*c\b', low):
            gene = 'hzsC'
    genus = next((g for g in GENERA if g.lower() in low), None)
    acc = h.split()[0]
    genus = genus or KNOWN_GENUS.get(acc)
    return gene, genus


def eutils(path, **params):
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/%s?%s' % (path, urllib.parse.urlencode(params))
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read().decode()
        except Exception as e:
            log('  E-utilities falhou (%s), tentando de novo' % e); time.sleep(2)
    return ''


def step_refs():
    refs = {}
    if os.path.exists(LOCAL_REFS):
        for h, s in read_fasta(LOCAL_REFS):
            gene, genus = classify_header(h)
            if gene and genus:
                refs[h.split()[0]] = (h, s, gene, genus, 'local')
    log('referencias locais anammox:', len(refs))
    for genus in GENERA:
        for gene, terms in GENE_TERMS.items():
            q = '(%s) AND "Candidatus %s"[Organism]' % (' OR '.join('"%s"[Title]' % t for t in terms), genus)
            j = eutils('esearch.fcgi', db='protein', term=q, retmax='4', retmode='text')
            ids = re.findall(r'<Id>(\d+)</Id>', j)
            time.sleep(0.4)
            if not ids:
                continue
            txt = eutils('efetch.fcgi', db='protein', id=','.join(ids), rettype='fasta', retmode='text')
            time.sleep(0.4)
            for blk in txt.split('>')[1:]:
                h, *sq = blk.split('\n')
                g2, gn2 = classify_header(h)
                if g2 == gene and gn2 == genus:
                    refs[h.split()[0]] = (h, ''.join(sq).strip(), g2, gn2, 'NCBI')
    tab = pd.DataFrame([(a, v[2], v[3], v[4], v[0]) for a, v in refs.items()],
                       columns=['accession', 'gene', 'genus', 'source', 'header'])
    tab = tab.sort_values(['gene', 'genus', 'accession'])
    tab.to_csv(f'{OUT}/anammox_reference_panel.tsv', sep='\t', index=False)
    write_fasta(f'{OUT}/anammox_refs.faa', [(a, v[1]) for a, v in refs.items()])
    log('painel de referencias: %d sequencias' % len(refs)); print(tab.groupby(['gene', 'genus']).size().unstack(fill_value=0).to_string())
    return tab


def step_ref_identity(orfs, tab):
    refs = f'{OUT}/anammox_refs.faa'
    run([need('diamond'), 'makedb', '--in', refs, '-d', f'{OUT}/anammox_refs', '--quiet'])
    tsv = f'{OUT}/orfs_vs_anammox_refs.tsv'
    run([need('diamond'), 'blastp', '-q', orfs, '-d', f'{OUT}/anammox_refs', '-o', tsv, '--more-sensitive', '-k', '200',
         '-e', '1e-10', '-p', str(THREADS), '--quiet', '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length',
         'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
    h = pd.read_csv(tsv, sep='\t', header=None, names=['orf', 'acc', 'pident', 'alnlen', 'qlen', 'slen', 'qs', 'qe', 'ss',
                                                        'se', 'evalue', 'bitscore'])
    h = h.merge(tab[['accession', 'gene', 'genus']], left_on='acc', right_on='accession')
    h['qcov'] = (h.qe - h.qs + 1) / h.qlen
    h['scov'] = (h.se - h.ss + 1) / h.slen
    h = h[(h.pident >= 40)]
    rows = []
    for orf, g in h.groupby('orf'):
        gene = g.sort_values('bitscore', ascending=False).iloc[0].gene
        gg = g[g.gene == gene]
        per = gg.sort_values('pident', ascending=False).groupby('genus').first()
        per = per.sort_values('pident', ascending=False)
        best = per.iloc[0]
        second = per.iloc[1] if len(per) > 1 else None
        rows.append(dict(orf=orf, contig=orf.rsplit('_', 1)[0], gene=gene, orf_len_aa=int(gg.qlen.iloc[0]),
                         best_genus=per.index[0], best_ref=best.acc, best_pident=round(best.pident, 1),
                         best_scov=round(best.scov, 2),
                         second_genus=per.index[1] if second is not None else '',
                         second_pident=round(second.pident, 1) if second is not None else '',
                         margin=round(best.pident - second.pident, 1) if second is not None else '',
                         all_genera=';'.join('%s:%.1f' % (k, v.pident) for k, v in per.iterrows())))
    df = pd.DataFrame(rows).sort_values(['gene', 'orf'])
    df.to_csv(f'{OUT}/orf_identity_by_genus.tsv', sep='\t', index=False)
    return df


# ---------------------------------------------------------------- 5. arvores
def step_trees(orfs, ident, tab):
    mafft, trimal = shutil.which('mafft'), shutil.which('trimal')
    iq = next((shutil.which(x) for x in ('iqtree3', 'iqtree2', 'iqtree') if shutil.which(x)), None)
    if not (mafft and trimal and iq):
        log('AVISO: mafft/trimal/iqtree nao estao todos no PATH; arvores puladas (veja o fim do script para instalar).')
        return
    orfseq = dict(read_fasta(orfs)); refseq = dict(read_fasta(f'{OUT}/anammox_refs.faa'))
    for gene in ['hzsA', 'hdh']:
        d = f'{OUT}/tree_{gene}'; os.makedirs(d, exist_ok=True)
        recs = []
        for _, r in ident[ident.gene == gene].iterrows():
            recs.append(('THISSTUDY_%s' % r.orf, orfseq[r.orf]))
        for _, r in tab[tab.gene == gene].iterrows():
            recs.append(('%s_%s' % (r.genus, r.accession.replace('.', 'v')), refseq[r.accession]))
        if len(recs) < 6:
            log('poucas sequencias para', gene); continue
        write_fasta(f'{d}/in.faa', recs)
        with open(f'{d}/aln.faa', 'w') as f:
            run([mafft, '--auto', '--thread', str(THREADS), f'{d}/in.faa'], stdout=f)
        run([trimal, '-in', f'{d}/aln.faa', '-out', f'{d}/aln_trim.faa', '-automated1'])
        run([iq, '-s', f'{d}/aln_trim.faa', '-m', 'MFP', '-B', '1000', '-alrt', '1000', '-T', 'AUTO', '--prefix',
             f'{d}/{gene}', '-redo', '-quiet'])
        log('arvore:', f'{d}/{gene}.treefile')


# ---------------------------------------------------------------- resumo
def summary(pg, ident):
    lines = []
    for c in CONTIGS:
        sub = pg[pg.contig == c]
        gen = sub.genus.replace('', '?').value_counts().to_dict()
        ph = sub.phylum.replace('', '?').value_counts().to_dict()
        lines.append('%s: %d ORFs com hit no proGenomes; generos do melhor hit: %s; filos: %s' % (c, len(sub), gen, ph))
        for _, r in ident[ident.contig == c].iterrows():
            lines.append('    %s (%s, %d aa): melhor genero %s %.1f%% (%s), segundo %s %s%%, margem %s' % (
                r.orf, r.gene, r.orf_len_aa, r.best_genus, r.best_pident, r.best_ref, r.second_genus, r.second_pident, r.margin))
    txt = '\n'.join(lines)
    open(f'{OUT}/SUMMARY.txt', 'w').write(txt + '\n')
    print('\n' + txt)


if __name__ == '__main__':
    fna = step_contigs()
    orfs = step_orfs(fna)
    pg = step_progenomes(orfs)
    tab = step_refs()
    ident = step_ref_identity(orfs, tab)
    step_trees(orfs, ident, tab)
    summary(pg, ident)
    log('FIM. Saidas em', OUT)
