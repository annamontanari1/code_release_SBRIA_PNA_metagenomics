# Regras de evidência compartilhadas pela Figure 1 (auditoria) e pela Figure S2.
# Critério (definido pela autora, 16/10/2026):
#   - hit considerado: e-value <= 1e-10; melhor hit = maior bit score
#   - supported: identidade >= 40% E cobertura >= 70% do comprimento real da
#     referência (fastas de results/targeted_search/) E, para famílias com
#     parálogos próximos, anotação DRAM/Pfam concordante no gene/contig do hit
#   - partial: passa identidade e anotação, mas cobertura < 70%
#   - genes partidos em ORFs adjacentes do mesmo contig: soma-se a cobertura se
#     os trechos da referência não se sobrepõem (coordenadas em
#     results/targeted_search/narG_bin008_ORF_refcoords.tsv, blastp das ORFs
#     existentes contra a referência existente) e o hit não é marcado parcial.
import csv, re, os, collections

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
TS = os.path.join(ROOT, 'results', 'targeted_search') + os.sep
EV_MAX, ID_MIN, COV_MIN = 1e-10, 40.0, 70.0

def _fasta_len(fn):
    d = {}; k = None
    for l in open(fn):
        l = l.strip()
        if l.startswith('>'): k = l[1:].split()[0]; d[k] = 0
        elif k: d[k] += len(l)
    return d

REFLEN = {}
for f in ['amo_refs.fasta', 'hzs_hdh_refs.fasta', 'nxr_refs.fasta',
          'denitrificacao_refs.fasta', 'nirS_ref.fasta', 'paracoccus_nirS.fasta',
          'nxr_anammox_ref.fasta']:  # Ca. K. stuttgartiensis NxrA (kustd1700), 16/10/2026
    REFLEN.update(_fasta_len(TS + f))

# DRAM (leitura apenas): gene -> (KO, descrição KEGG, Pfam); contig -> genes
ANN, NODE_GENES, C2BIN = {}, collections.defaultdict(list), {}
KOS = collections.defaultdict(set)
with open(os.path.join(ROOT, 'DRAM annotation', 'annotations.tsv')) as f:
    r = csv.reader(f, delimiter='\t'); next(r)
    for row in r:
        b = row[1].split('.fasta')[0]; node = row[2].split('_length')[0]
        ANN[row[0]] = (row[8], row[9], row[17])
        NODE_GENES[node].append(row[0]); C2BIN[node] = b
        if row[8]: KOS[b].add(row[8])

# Coordenadas na referência de ORFs adjacentes (apenas pares encontrados)
ORF_COORDS = collections.defaultdict(list)
fn = TS + 'narG_bin008_ORF_refcoords.tsv'
if os.path.exists(fn):
    for l in open(fn):
        p = l.rstrip('\n').split('\t')
        ORF_COORDS[(p[0], re.sub(r'_\d+$', '', p[1]))].append((int(p[4]), int(p[5]), p[1]))

def load_hits(fn):
    out = []
    for i, l in enumerate(open(TS + fn), 1):
        p = l.rstrip('\n').split('\t')
        subj = p[1].split(' ')[0]
        out.append(dict(ref=p[0], subj=subj, node=re.search(r'(NODE_\d+)', subj).group(1),
                        pid=float(p[2]), aln=int(p[3]), ev=float(p[4]), bits=float(p[5]),
                        file=fn, line=i, is_gene=subj in ANN))
    return out

# ---- famílias com parálogos próximos ------------------------------------------
# conc: KOs concordantes; pf_ok: Pfam aceito quando o gene não tem KO;
# para: KOs que indicam o parálogo (anotação discordante).
FAMILIES = {
    'hdh':  dict(conc={'K20935'}, pf_ok=set(), para={'K10535', 'K00378'}),
    'napA': dict(conc={'K02567'}, pf_ok=set(), para=None),       # qualquer outro KO = discordante
    'narG/nxrA family': dict(conc={'K00370'}, pf_ok={'PF00384'}, para=None),
    'nxrA': dict(conc={'K00370'}, pf_ok={'PF00384'}, para=None),
    'amoC': dict(conc={'K10946'}, pf_ok={'PF04896'}, para=None),
}

def _gene_ok(g, fam):
    ko, desc, pf = ANN[g]
    if ko: return ('yes' if ko in fam['conc'] else 'no'), f'{ko} {desc[:45]}'
    if any(x in pf for x in fam['pf_ok']): return 'yes', 'Pfam ' + ','.join(x for x in fam['pf_ok'] if x in pf)
    return 'none', ''

def family_check(gene, hit):
    """Retorna (aceito, nota). Sem anotação DRAM no locus = aceito com nota."""
    fam = FAMILIES.get(gene)
    if fam is None: return True, ''
    if hit['is_gene']:
        st, txt = _gene_ok(hit['subj'], fam)
        if st == 'no': return False, f'paralog: hit gene annotated {txt}'
        return True, (f'annotation concordant ({txt})' if st == 'yes' else 'hit gene without KO/Pfam')
    genes = NODE_GENES.get(hit['node'], [])
    if not genes: return True, 'unbinned contig, no DRAM annotation'
    res = [(g,) + _gene_ok(g, fam) for g in genes]
    conc = [r for r in res if r[1] == 'yes']
    if conc: return True, f'annotation concordant on contig ({conc[0][2]})'
    # hit no nível de contig (tBLASTn): o gene exato não é conhecido, então só
    # um KO de parálogo conhecido no contig conta como discordante
    paras = [r for r in res if fam['para'] and ANN[r[0]][0] in fam['para']]
    if paras: return False, f'paralog: contig gene annotated {paras[0][2]}'
    return True, 'no concordant or paralog annotation on contig'

def coverage(hit, all_hits):
    """Cobertura (%) da referência; soma ORFs adjacentes sem sobreposição."""
    base = 100 * hit['aln'] / REFLEN[hit['ref']]
    if not hit['is_gene']: return base, ''
    contig = re.sub(r'_\d+$', '', hit['subj']); n = int(hit['subj'].rsplit('_', 1)[1])
    adj = [h for h in all_hits if h['ref'] == hit['ref'] and h['is_gene'] and h['ev'] <= EV_MAX
           and re.sub(r'_\d+$', '', h['subj']) == contig and abs(int(h['subj'].rsplit('_', 1)[1]) - n) == 1]
    if not adj: return base, ''
    co = sorted(ORF_COORDS.get((hit['ref'], contig), []))
    if len(co) < 2: return base, 'adjacent ORF hit, no reference coordinates (not summed)'
    if any(b[0] <= a[1] for a, b in zip(co, co[1:])):
        return base, 'adjacent ORFs overlap on reference (not summed)'
    tot = sum(e - s + 1 for s, e, _ in co)
    orfs = '+'.join('_' + c[2].rsplit('_', 1)[1] for c in co)
    return 100 * tot / REFLEN[hit['ref']], f'split gene: ORFs {orfs} summed ({"; ".join(f"{s}-{e}" for s, e, _ in co)})'

def evaluate(gene, cands, all_hits):
    """cands: hits da unidade para a referência. Retorna dict com status e motivo."""
    cands = [h for h in cands if h['ev'] <= EV_MAX]
    if not cands: return dict(status='no hit', best=None, cov=None, note='')
    best = max(cands, key=lambda h: h['bits'])
    cov, cnote = coverage(best, all_hits)
    ok, fnote = family_check(gene, best)
    notes = [x for x in (cnote, fnote) if x]
    if best['pid'] < ID_MIN:
        return dict(status='not supported', best=best, cov=cov,
                    note='; '.join([f'identity {best["pid"]:.1f}% < 40%'] + notes))
    if not ok:
        return dict(status='not supported', best=best, cov=cov, note='; '.join(notes))
    return dict(status='supported' if cov >= COV_MIN else 'partial', best=best, cov=cov,
                note='; '.join(notes))
