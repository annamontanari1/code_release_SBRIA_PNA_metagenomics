"""Junta CheckM v1.0.18 (Table S1), CheckM2 e GUNC por bin; categorias MIMAG; cópia da Table S1 com colunas novas.
Rodar com o python do ambiente gunc (usa as funções do próprio GUNC para reconstruir a divisão de clados)."""
import os, csv, collections
import pandas as pd
from gunc.get_scores import read_diamond_output, create_base_data, get_abundant_lineages_cutoff

QC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.abspath(os.path.join(QC, '..', '..', '..'))
S1 = os.path.join(ROOT, 'results', 'tables_text', 'TableS1_bins.tsv')


def cat_cc(comp, cont):
    """MIMAG só por completude/contaminação."""
    if comp >= 90 and cont < 5: return 'High'
    if comp >= 50 and cont < 10: return 'Medium'
    return 'Low'


def cat_full(comp, cont, rrna, ntrna):
    """Mesma regra de make_TableS1.py: High exige 5S/16S/23S presentes e >= 18 tipos de tRNA."""
    allr = all(x.strip() not in ('–', '-', '') for x in rrna.split('/'))
    if comp >= 90 and cont < 5 and allr and ntrna >= 18: return 'High'
    return 'Medium' if comp >= 50 and cont < 10 else 'Low'


s1 = pd.read_csv(S1, sep='\t')
cm2 = pd.read_csv(os.path.join(QC, 'checkm2_out', 'quality_report.tsv'), sep='\t').set_index('Name')
gdir = os.path.join(QC, 'gunc_out')
gmax = pd.read_csv([os.path.join(gdir, f) for f in os.listdir(gdir) if f.endswith('maxCSS_level.tsv')][0], sep='\t').set_index('genome')

# divisão de clados no nível do CSS máximo, com o mesmo corte (> 2% dos genes mapeados) do GUNC
split = {}
dd = os.path.join(gdir, 'diamond_output')
for f in sorted(os.listdir(dd)):
    g = f.split('.diamond.')[0]
    df = read_diamond_output(os.path.join(dd, f))
    base = create_base_data(df, 'progenomes_2.1', None)
    lvl = gmax.loc[g, 'taxonomic_level']
    cut = get_abundant_lineages_cutoff(False, len(df))
    cnt = base[lvl].value_counts()
    cnt = cnt[cnt > cut]
    tot = cnt.sum()
    split[g] = '; '.join(f"{k.split(' ', 1)[-1]} {100 * v / tot:.1f}%" for k, v in cnt.items())

rows = []
for _, r in s1.iterrows():
    b = r['Bin']; c1, k1 = r['Completeness (%)'], r['Contamination (%)']
    c2, k2 = cm2.loc[b, 'Completeness'], cm2.loc[b, 'Contamination']
    G = gmax.loc[b]
    rr, nt = r['rRNA genes (5S/16S/23S)'], int(r['tRNA types'])
    m1, m2 = cat_cc(c1, k1), cat_cc(c2, k2)
    f1, f2 = cat_full(c1, k1, rr, nt), cat_full(c2, k2, rr, nt)
    rows.append(dict(bin=b, taxonomy=r['GTDB-Tk taxonomy (R214)'].split(';g__')[0].split(';')[-1] + ';g__' + r['GTDB-Tk taxonomy (R214)'].split(';g__')[-1].split(';')[0],
                     contigs=r['Contigs'], size_Mb=r['Genome size (Mb)'],
                     CheckM1_completeness=c1, CheckM1_contamination=k1,
                     CheckM2_completeness=c2, CheckM2_contamination=k2,
                     CheckM2_model=cm2.loc[b, 'Completeness_Model_Used'],
                     d_completeness=round(c2 - c1, 2), d_contamination=round(k2 - k1, 2),
                     GUNC_pass=G['pass.GUNC'], GUNC_CSS=round(G['clade_separation_score'], 3),
                     GUNC_contamination_portion=round(G['contamination_portion'], 3),
                     GUNC_n_effective_surplus_clades=round(G['n_effective_surplus_clades'], 2),
                     GUNC_RRS=round(G['reference_representation_score'], 3),
                     GUNC_level_maxCSS=G['taxonomic_level'], GUNC_clade_split_at_level=split.get(b, ''),
                     GUNC_genes_mapped=G['n_genes_mapped'], GUNC_genes_called=G['n_genes_called'],
                     rRNA=rr, tRNA_types=nt,
                     MIMAG_cc_CheckM1=m1, MIMAG_cc_CheckM2=m2, MIMAG_full_CheckM1=f1, MIMAG_full_CheckM2=f2,
                     category_changes=('yes: ' + (f'cc {m1}->{m2}' if m1 != m2 else '') + (' ' if m1 != m2 and f1 != f2 else '') + (f'full {f1}->{f2}' if f1 != f2 else '')) if (m1 != m2 or f1 != f2) else 'no'))
out = pd.DataFrame(rows)
out.to_csv(os.path.join(QC, 'QC_comparison_CheckM1_CheckM2_GUNC.tsv'), sep='\t', index=False)

# cópia da Table S1 com colunas novas (a original não é alterada)
s1c = s1.copy()
s1c.insert(s1c.columns.get_loc('Contamination (%)') + 1, 'CheckM2 completeness (%)', [cm2.loc[b, 'Completeness'] for b in s1c['Bin']])
s1c.insert(s1c.columns.get_loc('CheckM2 completeness (%)') + 1, 'CheckM2 contamination (%)', [cm2.loc[b, 'Contamination'] for b in s1c['Bin']])
s1c.insert(s1c.columns.get_loc('CheckM2 contamination (%)') + 1, 'GUNC (pass/fail)', ['pass' if gmax.loc[b, 'pass.GUNC'] in (True, 'True') else 'fail' for b in s1c['Bin']])
s1c.to_csv(os.path.join(QC, 'TableS1_bins_with_CheckM2_GUNC.tsv'), sep='\t', index=False)
print(out[['bin', 'CheckM1_completeness', 'CheckM1_contamination', 'CheckM2_completeness', 'CheckM2_contamination', 'GUNC_pass', 'GUNC_CSS',
           'GUNC_contamination_portion', 'GUNC_level_maxCSS', 'MIMAG_full_CheckM1', 'MIMAG_full_CheckM2', 'category_changes']].to_string())
