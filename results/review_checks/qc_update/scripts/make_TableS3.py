"""Table S3 (QC): CheckM v1.0.18 (Table S1) x CheckM2 x GUNC, uma linha por bin. Não altera a Table S1."""
import os
import pandas as pd
from decimal import Decimal, ROUND_HALF_UP
def r1(x, n=1):
    """Arredonda meio-para-cima sobre o valor decimal (94.85 -> 94.9), não sobre o binário."""
    return str(Decimal(str(x)).quantize(Decimal(1).scaleb(-n), rounding=ROUND_HALF_UP))
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment

QC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.abspath(os.path.join(QC, '..', '..', '..'))
s1 = pd.read_csv(os.path.join(ROOT, 'results', 'tables_text', 'TableS1_bins.tsv'), sep='\t').set_index('Bin')
cm2 = pd.read_csv(os.path.join(QC, 'checkm2_out', 'quality_report.tsv'), sep='\t').set_index('Name')
g = pd.read_csv(os.path.join(QC, 'gunc_out', 'GUNC.progenomes_2.1.maxCSS_level.tsv'), sep='\t').set_index('genome')

RANK = {'g': 'genus', 'f': 'family', 'o': 'order', 'c': 'class', 'p': 'phylum', 'd': 'domain'}
def lowest(tax):
    """Gênero quando resolvido; senão o menor nível resolvido, com o nível entre parênteses."""
    for r in ['g', 'f', 'o', 'c', 'p', 'd']:
        for t in tax.split(';'):
            if t.startswith(r + '__') and t[3:]:
                return t[3:] if r == 'g' else f'{t[3:]} ({RANK[r]})'

NOTE = {
    'bin.001': 'CheckM2 completeness below 50%; MIMAG category changes from medium to low quality.',
    'bin.002': 'Poorly represented in the GUNC reference (RRS 0.25); GUNC scores should be interpreted with caution.',
    'bin.004': 'Fails GUNC at phylum level: 8 contigs (81 kb) assigned to Proteobacteria, at lower depth than the rest of the bin (53× vs 66×).',
    'bin.005': 'Poorly represented in the GUNC reference (RRS 0.26); GUNC scores should be interpreted with caution.',
    'bin.006': 'Fails GUNC at genus level (6.5% of genes assigned to Ca. Accumulibacter).',
    'bin.007': 'CheckM2 contamination below 10%; MIMAG category changes from low to medium quality.',
    'bin.012': 'Higher contamination with CheckM2; passes GUNC, with most non-Armatimonadota genes on two contigs (NODE_165, NODE_508).',
    'bin.015': 'Chimeric by CheckM, CheckM2 and GUNC (contigs assigned to Planctomycetes, 10.7% of genes).',
    'bin.018': 'GTDB-Tk multiple-marker warning not supported by GUNC (low CSS).',
    'bin.023': 'Poorly represented in the GUNC reference (RRS 0.23); GUNC scores should be interpreted with caution.',
}
HDR = ['Bin', 'GTDB-Tk (genus or lowest resolved rank)', 'CheckM v1.0.18 completeness (%)', 'CheckM v1.0.18 contamination (%)',
       'CheckM2 completeness (%)', 'CheckM2 contamination (%)', 'GUNC', 'GUNC CSS', 'GUNC contamination portion', 'Note']
rows = []
for b in [f'bin.{i:03d}' for i in range(1, 26)]:
    rows.append([b, lowest(s1.loc[b, 'GTDB-Tk taxonomy (R214)']),
                 r1(s1.loc[b, 'Completeness (%)']), r1(s1.loc[b, 'Contamination (%)']),
                 r1(cm2.loc[b, 'Completeness']), r1(cm2.loc[b, 'Contamination']),
                 'pass' if g.loc[b, 'pass.GUNC'] in (True, 'True') else 'fail',
                 r1(g.loc[b, 'clade_separation_score'], 2), r1(g.loc[b, 'contamination_portion'], 2), NOTE.get(b, '')])
with open(os.path.join(QC, 'TableS3_QC.tsv'), 'w') as w:
    w.write('\t'.join(HDR) + '\n')
    for r in rows: w.write('\t'.join(r) + '\n')

wb = Workbook(); ws = wb.active; ws.title = 'Table S3'
ws.append(HDR)
for c in ws[1]: c.font = Font(bold=True); c.alignment = Alignment(wrap_text=True, vertical='top')
for r in rows:
    # números como número no Excel, com formato de 1 ou 2 casas
    ws.append([r[0], r[1]] + [float(x) for x in r[2:6]] + [r[6], float(r[7]), float(r[8]), r[9]])
for row in ws.iter_rows(min_row=2):
    for c in row[2:6]: c.number_format = '0.0'
    for c in row[7:9]: c.number_format = '0.00'
    row[9].alignment = Alignment(wrap_text=True, vertical='top')
for col, wdt in zip('ABCDEFGHIJ', [9, 30, 14, 14, 13, 13, 8, 9, 13, 60]): ws.column_dimensions[col].width = wdt
wb.save(os.path.join(QC, 'TableS3_QC.xlsx'))
print(pd.read_csv(os.path.join(QC, 'TableS3_QC.tsv'), sep='\t', dtype=str).iloc[:, :9].to_string(index=False))
