import os, re, csv, collections
ROOT=os.path.expanduser('~/sbria_anammox_shotgun'); RC=ROOT+'/results/review_checks'
D=ROOT+'/das_tool_result (1)'
def node(x): return re.match(r'(NODE_\d+)',x).group(1)
# contig info
INFO={}
for l in open(RC+'/00_base/contig_len_gc.tsv'):
    n,L,g=l.rstrip('\n').split('\t'); INFO[node(n)]=dict(full=n,len=int(L),gc=float(g),spades_cov=float(n.split('_cov_')[1]))
for r in csv.reader(open(ROOT+'/results/contig_check/asm_coverage.tsv'),delimiter='\t'):
    if r[0].startswith('#'): continue
    k=node(r[0]); INFO[k]['depth']=float(r[6]); INFO[k]['numreads']=int(r[3]); INFO[k]['mapq']=float(r[8])
KEY={}
for r in csv.reader(open(D+'/das_tool_output_dir_DASTool_bins/das_tool_name_key.tsv'),delimiter='\t'):
    if r[0].startswith('Original'): continue
    KEY[r[0].replace('.fa','')]=r[1].replace('.fasta','')
FINAL={}
for l in open(D+'/das_tool_output_dir_DASTool_scaffolds2bin.txt'):
    c,b=l.split(); FINAL[node(c)]=KEY[b]+' ('+b.replace('NGS1140_6S_','')+')'
FINALBIN={k:v.split(' ')[0] for k,v in FINAL.items()}
BINNER=collections.defaultdict(dict)
for t in ['metabat2','maxbin2','concoct']:
    for l in open(f'{D}/NGS1140_6S_{t}_bins_contigs_to_bins.tsv'):
        c,b=l.split()[:2]; BINNER[node(c)][t]=b.replace(f'NGS1140_6S_{t}_bins.','')
def binners(n): 
    d=BINNER.get(n,{}); return '; '.join(f'{t}:{d[t]}' for t in ['metabat2','maxbin2','concoct'] if t in d) or 'none'
REFLEN={l.split()[0]:int(l.split()[1]) for l in open(RC+'/refs/refs_len.tsv')}
REFDESC={}
for l in open(RC+'/refs/refs_clean.faa'):
    if l[0]=='>': a=l[1:].strip().split(' ',1); REFDESC[a[0]]=a[1] if len(a)>1 else ''
def bin_contigs(b):
    return [k for k,v in FINALBIN.items() if v==b]
