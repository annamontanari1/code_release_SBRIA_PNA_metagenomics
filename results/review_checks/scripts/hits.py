import sys; sys.path.insert(0,__file__.rsplit('/',1)[0]); from base import *
def load_tbn(ev=1e-10):
    H=[]
    for l in open(RC+'/00_base/tblastn_allrefs_vs_assembly.tsv'):
        p=l.rstrip('\n').split('\t')
        if float(p[10])>ev: continue
        H.append(dict(ref=p[0],contig=node(p[1]),pid=float(p[2]),aln=int(p[3]),qs=int(p[6]),qe=int(p[7]),ss=int(p[8]),se=int(p[9]),ev=float(p[10]),bits=float(p[11]),qlen=int(p[12])))
    return H
def load_bp(ev=1e-10):
    H=[]
    for l in open(RC+'/00_base/blastp_allrefs_vs_hitORFs.tsv'):
        p=l.rstrip('\n').split('\t')
        if float(p[10])>ev: continue
        H.append(dict(ref=p[0],orf=p[1],contig=node(p[1]),pid=float(p[2]),aln=int(p[3]),qs=int(p[6]),qe=int(p[7]),ss=int(p[8]),se=int(p[9]),ev=float(p[10]),bits=float(p[11]),qlen=int(p[12]),slen=int(p[13]),qcov=float(p[14]),scov=round(100*(int(p[9])-int(p[8])+1)/int(p[13]),1)))
    return H
def union(iv):
    iv=sorted(iv); tot=0; cs=ce=None
    for s,e in iv:
        if cs is None or s>ce+1:
            if cs is not None: tot+=ce-cs+1
            cs,ce=s,e
        else: ce=max(ce,e)
    if cs is not None: tot+=ce-cs+1
    return tot
def per_contig(H, ref):
    """best HSP per contig plus union of query coverage over all HSPs on that contig"""
    d=collections.defaultdict(list)
    for h in H:
        if h['ref']==ref: d[h['contig']].append(h)
    out=[]
    for c,hs in d.items():
        b=max(hs,key=lambda h:h['bits'])
        out.append(dict(b, n_hsp=len(hs), qcov_best=100*(b['qe']-b['qs']+1)/b['qlen'], qcov_union=100*union([(h['qs'],h['qe']) for h in hs])/b['qlen']))
    return sorted(out,key=lambda x:-x['bits'])
# ORF coordinates from prodigal headers
ORF={}
for l in open(RC+'/00_base/hit_contigs_orfs.faa'):
    if l[0]=='>':
        p=l[1:].split(' # '); ORF[p[0]]=(int(p[1]),int(p[2]),int(p[3]))
HDR='gene\tref\tref_len\tcontig\tcontig_len\tdepth\tGC\tfinal_bin\tbinners\tn_hsp\tpident_best\taln_len\tqstart\tqend\tsstart\tsend\tevalue\tbits\tqcov_best_pct\tqcov_union_pct\tORF\tORF_coords\tORF_len_aa\tORF_pid\tORF_qcov_pct\tORF_target_cov_pct\tORF_evalue\tORF_bits\tstatus\n'
def rows(gene, ref, H, BP, top=None, binfilter=None):
    out=[]
    for h in per_contig(H,ref)[:top]:
        c=h['contig']; I=INFO[c]; fb=FINALBIN.get(c,'unbinned')
        if binfilter and fb not in binfilter: continue
        bp=[b for b in BP if b['ref']==ref and b['contig']==c]; b=max(bp,key=lambda x:x['bits']) if bp else None
        cov=max(h['qcov_union'], b['qcov'] if b else 0)
        st='supported' if (h['pid']>=40 and cov>=70) else ('partial' if h['pid']>=40 else 'not supported (id<40%)')
        oc='%d-%d(%s)'%(ORF[b['orf']][0],ORF[b['orf']][1],'+' if ORF[b['orf']][2]==1 else '-') if b else ''
        out.append([gene,ref,REFLEN[ref],c,I['len'],round(I['depth'],2),I['gc'],fb,binners(c),h['n_hsp'],h['pid'],h['aln'],h['qs'],h['qe'],h['ss'],h['se'],h['ev'],h['bits'],round(h['qcov_best'],1),round(h['qcov_union'],1),
            b['orf'].split('_length')[0]+'_'+b['orf'].rsplit('_',1)[1] if b else '',oc,b['slen'] if b else '',b['pid'] if b else '',b['qcov'] if b else '',b['scov'] if b else '',b['ev'] if b else '',b['bits'] if b else '',st])
    return out
def write(fn, R):
    with open(fn,'w') as w:
        w.write(HDR)
        for r in R: w.write('\t'.join(map(str,r))+'\n')
