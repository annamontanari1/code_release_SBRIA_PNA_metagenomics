import csv,re,collections,glob
from openpyxl import Workbook
# 1. completeness/contamination from Fig1 script (CheckM v1.0.18)
cm={}
for l in open('Fig1_MAG_quality.R'):
    m=re.match(r'\s*"(bin\.\d+)",\s*"([^"]*)",\s*"([^"]*)",\s*([\d.]+),\s*([\d.]+)',l)
    if m: cm[m.group(1)]=(float(m.group(4)),float(m.group(5)))
# 2. GTDB-Tk taxonomy
tax={}
for fn in ['bacteria.csv','archea.csv']:
    for r in csv.DictReader(open(fn)):
        tax[r['User Genome'].replace('.fa','')]=r['Classification']
# 3. contig lengths and samtools depth (Bowtie2 + samtools coverage)
cov={}
for i,l in enumerate(open('results/contig_check/asm_coverage.tsv')):
    if i==0: continue
    p=l.rstrip('\n').split('\t'); cov[p[0].split()[0]]=(int(p[2]),float(p[6]))
# 4. DAS Tool scaffolds2bin (all contigs per DAS bin)
s2b=collections.defaultdict(list)
for l in open('das_tool_result (1)/das_tool_output_dir_DASTool_scaffolds2bin.txt'):
    p=l.split()
    if len(p)>=2: s2b[p[1]].append(p[0])
# 5. DRAM bin -> contig set (>=2500 bp) -> dominant DAS bin
dram=collections.defaultdict(set)
with open('DRAM annotation/annotations.tsv') as f:
    for row in csv.DictReader(f,delimiter='\t'):
        m=re.match(r'(NODE_\d+_length_\d+_cov_[\d.]+)',row['scaffold'])
        if m: dram[row['fasta'].split('.fasta')[0]].add(m.group(1))
c2das={c:b for b,cs in s2b.items() for c in cs}
# 6. rRNA (barrnap) e tRNA (tRNAscan-SE) do DRAM -> MIMAG completo (16/10/2026)
part=set()
for r in csv.DictReader(open('DRAM annotation/rrnas.tsv'),delimiter='\t'):
    if r['note'].strip(): part.add((r['fasta'].split('.fasta')[0],r['type'].split()[0]))
gs={r['genome'].split('.fasta')[0]:r for r in csv.DictReader(open('DRAM annotation/genome_stats.tsv'),delimiter='\t')}
trna=collections.defaultdict(set)
for r in csv.DictReader(open('DRAM annotation/trnas.tsv'),delimiter='\t'):
    ty=r['Type'].strip()
    if ty and ty not in ('Undet','Pseudo','SeC','Sup','fMet','Ile2'): trna[r['fasta'].split('.fasta')[0]].add(ty)
def rrna(b):
    out=[];allr=True
    for k in ('5S','16S','23S'):
        if gs[b][k+' rRNA'].strip(): out.append(k+(' p' if (b,k) in part else ''))
        else: out.append('–'); allr=False
    return '/'.join(out),allr
rows=[]
for b in sorted(cm):
    cnt=collections.Counter(c2das.get(c,'?') for c in dram[b])
    das,n=cnt.most_common(1)[0]
    frac=n/max(1,len(dram[b]))
    cs=[c for c in s2b[das] if c in cov]
    size=sum(cov[c][0] for c in cs)
    depth=sum(cov[c][0]*cov[c][1] for c in cs)/size
    comp,con=cm[b]
    rr,allr=rrna(b); nt=len(trna[b])
    tier='High' if comp>=90 and con<5 and allr and nt>=18 else ('Medium' if comp>=50 and con<10 else 'Low')
    t=tax.get(b+'.fa','') or tax.get(b,'')
    parts=[x for x in t.split(';')]
    ph=[x for x in parts if x.startswith('p__')]; g=[x for x in parts if x.startswith('g__')]; f_=[x for x in parts if x.startswith('f__')]
    s=[x for x in parts if x.startswith('s__')]
    rows.append([b,t,comp,con,round(size/1e6,2),len(cs),round(depth,1),rr,nt,tier,das,round(frac,2)])
hdr=['Bin','GTDB-Tk taxonomy (R214)','Completeness (%)','Contamination (%)','Genome size (Mb)','Contigs','Mean depth (x)','rRNA genes (5S/16S/23S)','tRNA types','Tier (MIMAG)','DAS Tool bin (internal)','Fraction of DRAM contigs matched']
with open('results/tables_text/TableS1_bins.tsv','w') as o:
    o.write('\t'.join(hdr)+'\n')
    for r in rows: o.write('\t'.join(map(str,r))+'\n')
wb=Workbook(); ws=wb.active; ws.title='Table S1'; ws.append(hdr)
for r in rows: ws.append(r)
wb.save('results/tables_text/TableS1_bins.xlsx')
for r in rows: print(r[0],r[2],r[3],r[6],r[7],r[8],r[9])
