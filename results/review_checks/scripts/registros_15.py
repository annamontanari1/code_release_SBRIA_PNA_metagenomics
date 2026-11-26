#!/usr/bin/env python3
"""Secao 15 do README (21/10/2026): verificacoes de registro para o manuscrito. Nao altera resultados existentes.

  1. trimAl/IQ-TREE da arvore de AmoA: versoes registradas + trimAl -automated1 refeito sobre o MAFFT original (copia em 15_registros/)
  2. pares de reads NODE_13511 -> outros contigs: mesmo criterio de link_hzs_contigs_to_bins.py / comammox/item1_contigs.py
     (leitura primaria e mapeada, mate mapeado em OUTRO contig; coluna conservadora = MAPQ >= 20 da leitura)
  3. sinal Tat de NODE_4298_3 / NODE_1742_7: so le os registros existentes (nenhuma ferramenta nova)
  4. Kaiju: metadados locais (relatorios do kb_kaiju)

Saidas: results/review_checks/15_registros/
Ambiente: ~/miniconda3/envs/phylo/bin/python (pysam); trimAl de ~/miniconda3/envs/trimal (o mesmo usado em 19/10/2026).
"""
import os, re, sys, subprocess, collections, datetime, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from base import *   # ROOT, RC, INFO, FINALBIN, node
import pysam

OUT = RC + '/15_registros/'; os.makedirs(OUT, exist_ok=True)
TREE = RC + '/comammox/amoA_tree/'
TRIMAL = os.path.expanduser('~/miniconda3/envs/trimal/bin/trimal')
MAPQ_MIN = 20
rep = []
def say(s=''): rep.append(s); print(s)
def md5(p): return hashlib.md5(open(p, 'rb').read()).hexdigest()
def mtime(p): return datetime.datetime.fromtimestamp(os.path.getmtime(p)).strftime('%Y-%m-%d %H:%M:%S')
def ncols(fa):
    seqs, k = collections.defaultdict(str), None
    for l in open(fa):
        if l[0] == '>': k = l[1:].strip()
        else: seqs[k] += l.strip()
    return len(seqs), sorted({len(s) for s in seqs.values()})

# ---------- 1. trimAl / IQ-TREE
say('== 1. trimAl e IQ-TREE (arvore de AmoA)')
say('IQ-TREE (cabecalho de amoA_tree.iqtree): ' + open(TREE + 'amoA_tree.iqtree').readline().strip())
cmd = [l for l in open(TREE + 'amoA_tree.log') if l.startswith('Command:')]
say('IQ-TREE (amoA_tree.log) ' + (cmd[0].strip() if cmd else 'Command: (nao registrado)'))
say('trimAl: nenhuma saida da ferramenta no diretorio; versao so citada em comammox/README.md e amoA_tree_plot.R ("v1.5.rev1")')
ver = subprocess.run([TRIMAL, '--version'], capture_output=True, text=True).stdout.strip()
meta = [f for f in os.listdir(os.path.expanduser('~/miniconda3/envs/trimal/conda-meta')) if f.startswith('trimal-')]
say(f'trimAl do ambiente "trimal": {ver} (pacote {meta[0][:-5] if meta else "?"}); conda-meta {mtime(os.path.expanduser("~/miniconda3/envs/trimal/conda-meta/" + meta[0]))}')
say(f'amoA_aln.faa {mtime(TREE + "amoA_aln.faa")}; amoA_aln_trim.faa {mtime(TREE + "amoA_aln_trim.faa")}')
rerun = OUT + 'amoA_aln_trim_rerun.faa'
subprocess.run([TRIMAL, '-in', TREE + 'amoA_aln.faa', '-out', rerun, '-automated1'], check=True)
n0, c0 = ncols(TREE + 'amoA_aln.faa'); n1, c1 = ncols(TREE + 'amoA_aln_trim.faa'); n2_, c2 = ncols(rerun)
say(f'colunas: MAFFT {c0} ({n0} seq) -> trim original {c1} ({n1} seq); refeito {c2} ({n2_} seq)')
say('arquivo refeito ' + ('IDENTICO byte a byte ao original' if md5(rerun) == md5(TREE + 'amoA_aln_trim.faa') else 'DIFERENTE do original')
    + f' (md5 {md5(rerun)})')
say()

# ---------- 2. pares de reads de NODE_13511
say('== 2. Pares de reads NODE_13511 -> outros contigs (BAM results/contig_check/asm.sorted.bam)')
bam = pysam.AlignmentFile(ROOT + '/results/contig_check/asm.sorted.bam')
n2c = {r: node(r) for r in bam.references}
C = 'NODE_13511'
tot = unm = same = 0
L, L20, Lboth, names = collections.Counter(), collections.Counter(), collections.Counter(), collections.defaultdict(set)
for r in bam.fetch(INFO[C]['full']):
    if r.is_secondary or r.is_supplementary or r.is_unmapped: continue
    tot += 1
    if r.mate_is_unmapped: unm += 1; continue
    m = n2c[r.next_reference_name]
    if m == C: same += 1; continue
    L[m] += 1; names[m].add(r.query_name)
    if r.mapping_quality >= MAPQ_MIN:
        L20[m] += 1
        if r.has_tag('MQ') and r.get_tag('MQ') >= MAPQ_MIN: Lboth[m] += 1
has_mq = any(True for r in bam.fetch(INFO[C]['full']) if r.has_tag('MQ'))
say(f'leituras primarias mapeadas: {tot}; mate nao mapeado: {unm}; mate no mesmo contig: {same}; mate em outro contig: {sum(L.values())}')
say(f'pares unicos (query_name) inter-contig: {sum(len(v) for v in names.values())}  (= leituras: cada par tem so uma leitura em {C})')
rows = []
for m, v in L.most_common():
    b = FINALBIN.get(m, 'unbinned')
    rows.append((C, m, b, INFO[m]['len'], round(INFO[m]['depth'], 2), v, L20[m], Lboth[m] if has_mq else 'NA'))
with open(OUT + 'NODE_13511_readpair_links.tsv', 'w') as w:
    w.write('contig\tmate_contig\tmate_final_bin\tmate_len\tmate_depth\tlinks_all\tlinks_read_MAPQ>=20\tlinks_both_MAPQ>=20\n')
    for x in rows: w.write('\t'.join(map(str, x)) + '\n')
for x in rows: say(f'  {x[1]} [{x[2]}; {x[3]} pb; {x[4]}x]: todos {x[5]}; MAPQ>=20 {x[6]}; ambos MAPQ>=20 {x[7]}')
byb, byb20 = collections.Counter(), collections.Counter()
for x in rows: byb[x[2]] += x[5]; byb20[x[2]] += x[6]
say('por bin (todos; MAPQ>=20): ' + '; '.join(f'{b}: {byb[b]} ({byb20[b]})' for b in byb))
say(f'total inter-contig: todos {sum(L.values())}; MAPQ>=20 {sum(L20.values())}')
say()

# ---------- 3. sinal Tat
say('== 3. Sinal Tat (NODE_4298_3, NODE_1742_7)')
for l in open(RC + '/03_mcrA_nxr/tat_Nterm50_motif.tsv'):
    if l.startswith(('NODE_4298', 'NODE_1742', 'Q1PZD8')): say('  ' + l.rstrip('\n').replace('\t', ' | '))
say('ferramenta: nenhuma (busca do motivo RR por regex nos 50 aa N-terminais; SignalP/TatP/TatFind/PRED-TAT/DeepTMHMM nao instalados nem citados em logs)')
say('origem da nota: results/targeted_search/nxr_bin009_vs_kuenenia_SUMMARY_not_used_by_scripts.tsv (16/10/2026), "N-terminal Tat signal (RR motif)"')
say()

# ---------- 4. Kaiju
say('== 4. Kaiju (banco RefSeq Genomes)')
K = ROOT + '/Kaiju/'
rp = K + 'kaiju_report/NGS1140_6S_trimmed_reads-genus.kaijuReport'
p = open(rp).readlines()[1].split('\t')[0]
job = re.search(r'output_(\d{13})', p)
say(f'arquivos: {", ".join(sorted(os.listdir(K)))}')
say(f'caminho registrado nos relatorios: {p}')
if job: say(f'carimbo do job KBase (epoch ms {job.group(1)}) = {datetime.datetime.fromtimestamp(int(job.group(1)) / 1000, datetime.timezone.utc):%Y-%m-%d %H:%M UTC}')
say(f'.kaiju: {mtime(K + "NGS1140_6S_trimmed_reads.kaiju")} (data do arquivo local)')
say('nodes.dmp / names.dmp / .fmi / log do kb_kaiju: nao encontrados localmente -> sem registro local da versao ou data do banco')
open(OUT + 'SUMMARY.txt', 'w').write('\n'.join(rep) + '\n')
