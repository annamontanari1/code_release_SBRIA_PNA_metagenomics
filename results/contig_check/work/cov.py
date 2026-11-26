import sys,re,collections
# lê SAM do stdin; cobertura = bases alinhadas (M/=/X) de alinhamentos primários com identidade >= 97% (NM/aln_len)
L={};cov=collections.Counter();n_all=collections.Counter();n_pass=collections.Counter()
cig=re.compile(r'(\d+)([MIDNSHP=X])')
tot=0
for line in sys.stdin:
    if line[0]=='@':
        if line.startswith('@SQ'):
            f=line.rstrip().split('\t');L[f[1][3:]]=int(f[2][3:])
        continue
    f=line.split('\t',11);flag=int(f[1])
    if flag&4 or flag&256 or flag&2048: continue
    tot+=1
    r=f[2];n_all[r]+=1
    nm=None
    if len(f)>11:
        for t in f[11].rstrip().split('\t'):
            if t.startswith('NM:i:'): nm=int(t[5:])
    al=sum(int(n) for n,o in cig.findall(f[5]) if o in 'MX=')
    if nm is None or al<50 or nm/al>0.03: continue
    n_pass[r]+=1;cov[r]+=al
print('contig\tlength\treads_aligned\treads_ge97id\tmean_depth_ge97id')
for r in sorted(L,key=lambda x:-L[x]):
    print('%s\t%d\t%d\t%d\t%.2f'%(r,L[r],n_all[r],n_pass[r],cov[r]/L[r]))
sys.stderr.write('total aligned primary records %d\n'%tot)
