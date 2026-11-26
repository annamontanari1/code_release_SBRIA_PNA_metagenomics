"""Contigs que carregam os clados minoritários no nível do CSS máximo (bins que falham no GUNC + destaques)."""
import os, re, sys
import pandas as pd
from gunc.get_scores import read_diamond_output, create_base_data, get_abundant_lineages_cutoff
QC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
cov = pd.read_csv(os.path.join(QC, '..', '..', 'contig_check', 'asm_coverage.tsv'), sep='\t')
cov['node'] = cov['#rname'].str.extract(r'(NODE_\d+)'); cov = cov.set_index('node')
gc = pd.read_csv(os.path.join(QC, '..', '00_base', 'contig_len_gc.tsv'), sep='\t', header=None, names=['c', 'len', 'gc'])
gc['node'] = gc['c'].str.extract(r'(NODE_\d+)'); gc = gc.set_index('node')
gmax = pd.read_csv(os.path.join(QC, 'gunc_out', 'GUNC.progenomes_2.1.maxCSS_level.tsv'), sep='\t').set_index('genome')
rows = []
for b in sys.argv[1:]:
    df = read_diamond_output(os.path.join(QC, 'gunc_out', 'diamond_output', f'{b}.diamond.progenomes_2.1.out'))
    base = create_base_data(df, 'progenomes_2.1', None)
    lvl = gmax.loc[b, 'taxonomic_level']; cut = get_abundant_lineages_cutoff(False, len(df))
    keep = base[lvl].value_counts(); keep = keep[keep > cut].index
    base = base[base[lvl].isin(keep)]
    major = base[lvl].value_counts().idxmax()
    for c, g in base.groupby('contig'):
        vc = g[lvl].value_counts(); n = node = re.search(r'(NODE_\d+)', c).group(1)
        top = vc.idxmax()
        rows.append(dict(bin=b, level=lvl, contig=n, contig_len=int(gc.loc[n, 'len']), depth=round(cov.loc[n, 'meandepth'], 2), GC=gc.loc[n, 'gc'],
                         genes_retained=int(vc.sum()), contig_majority_clade=top.split(' ', 1)[-1], frac_majority=round(vc.max() / vc.sum(), 2),
                         is_minor_clade_contig=top != major))
out = pd.DataFrame(rows)
out.to_csv(os.path.join(QC, 'GUNC_contig_clade_assignment_selected_bins.tsv'), sep='\t', index=False)
for b, g in out.groupby('bin'):
    m = g[g.is_minor_clade_contig]; M = g[~g.is_minor_clade_contig]
    w = lambda x, col: (x[col] * x.contig_len).sum() / x.contig_len.sum()
    print(f'{b}: {len(m)} minor-clade contigs ({m.contig_len.sum()} bp; {m.genes_retained.sum()} genes) depth {w(m,"depth") if len(m) else 0:.1f}x GC {w(m,"GC") if len(m) else 0:.1f} | major contigs depth {w(M,"depth"):.1f}x GC {w(M,"GC"):.1f}')
    print(m.sort_values('genes_retained', ascending=False).head(8)[['contig', 'contig_len', 'depth', 'GC', 'genes_retained', 'contig_majority_clade', 'frac_majority']].to_string(index=False))
