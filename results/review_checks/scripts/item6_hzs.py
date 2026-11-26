import sys; sys.path.insert(0,__file__.rsplit('/',1)[0]); from hits import *
import numpy as np, itertools, random, pysam, subprocess
random.seed(1)
OUT=RC+'/06_tableS2/'
HZS=['NODE_5668','NODE_5312','NODE_7170','NODE_17695','NODE_29475','NODE_23806','NODE_20914','NODE_47442','NODE_24842','NODE_24575']
# --- sequences
SEQ={}
def load(fn,keep):
    k=None
    for l in open(fn):
        if l[0]=='>':
            n=node(l[1:]); k=n if n in keep else None
            if k: SEQ[k]=[]
        elif k: SEQ[k].append(l.strip())
B9=bin_contigs('bin.009')
ALLB={b:bin_contigs(b) for b in sorted(set(FINALBIN.values()))}
need=set(HZS)|set(c for v in ALLB.values() for c in v)
load(ROOT+'/results/targeted_search/bin_assemblies/NGS1140_6S_assembly.fa',need)
SEQ={k:''.join(v).upper() for k,v in SEQ.items()}
# --- TNF (136 canonical tetramers)
comp=str.maketrans('ACGT','TGCA')
km=[''.join(p) for p in itertools.product('ACGT',repeat=4)]
can=sorted(set(min(k,k.translate(comp)[::-1]) for k in km)); idx={k:can.index(min(k,k.translate(comp)[::-1])) for k in km}
def tnf(s):
    v=np.zeros(len(can))
    for i in range(len(s)-3):
        j=idx.get(s[i:i+4]); 
        if j is not None: v[j]+=1
    return v/v.sum() if v.sum() else v
def centroid(cs):
    v=sum(tnf(SEQ[c])*len(SEQ[c]) for c in cs); return v/v.sum()
cent={b:centroid(cs) for b,cs in ALLB.items()}
b9cat=[SEQ[c] for c in B9]; b9w=[len(s) for s in b9cat]
def null(L,n=500):
    d=[]
    for _ in range(n):
        s=random.choices(b9cat,weights=b9w)[0]
        if len(s)<L: continue
        i=random.randrange(len(s)-L+1); d.append(np.linalg.norm(tnf(s[i:i+L])-cent['bin.009']))
    return np.array(d)
gcs=np.array([INFO[c]['gc'] for c in B9]); Ls=np.array([INFO[c]['len'] for c in B9])
gm=np.average(gcs,weights=Ls); gsd_w=np.sqrt(np.average((gcs-gm)**2,weights=Ls)); gsd_u=gcs.std(ddof=1); gmu=gcs.mean()
with open(OUT+'hzs_hdh_GC_TNF_vs_bin009.tsv','w') as w:
    w.write(f'# bin.009 GC: length-weighted {gm:.2f} ± {gsd_w:.2f} (SD); unweighted {gmu:.2f} ± {gsd_u:.2f} (n={len(B9)} contigs)\n')
    w.write('contig\tlen\tdepth\tGC\tGC_z_vs_bin009(weighted_sd)\tTNF_dist_to_bin009\tTNF_null_P95_same_len_windows_from_bin009\tTNF_percentile_vs_null\tnearest_bin_by_TNF\tTNF_dist_nearest\tTNF_rank_of_bin009_among_25\n')
    for c in HZS:
        v=tnf(SEQ[c]); d=np.linalg.norm(v-cent['bin.009']); nd=null(len(SEQ[c]))
        ds=sorted((np.linalg.norm(v-cent[b]),b) for b in cent); rank=[b for _,b in ds].index('bin.009')+1
        w.write('\t'.join(map(str,[c,INFO[c]['len'],round(INFO[c]['depth'],2),INFO[c]['gc'],round((INFO[c]['gc']-gm)/gsd_w,2),round(d,4),round(np.percentile(nd,95),4),round(100*(nd<d).mean(),1),ds[0][1],round(ds[0][0],4),rank]))+'\n')
# --- read-pair linkage
bam=pysam.AlignmentFile(ROOT+'/results/contig_check/asm.sorted.bam')
name2node={r:node(r) for r in bam.references}
with open(OUT+'hzs_hdh_readpair_links.tsv','w') as w, open(OUT+'hzs_hdh_readpair_links_detail.tsv','w') as wd:
    w.write('contig\tlen\tprimary_reads\treads_mate_unmapped\treads_mate_same_contig\treads_mate_other_contig\tlinks_to_bin009_contigs(all MAPQ)\tlinks_to_bin009_MAPQ>=20\tlinks_to_other_hzs_contigs\ttop_mate_contigs(n;bin)\n')
    wd.write('contig\tmate_contig\tmate_final_bin\tn_links_allMAPQ\tn_links_MAPQ>=20\n')
    for c in HZS:
        ref=INFO[c]['full']; n=un=same=oth=0; L=collections.Counter(); L20=collections.Counter()
        for r in bam.fetch(ref):
            if r.is_secondary or r.is_supplementary or r.is_unmapped: continue
            n+=1
            if r.mate_is_unmapped: un+=1; continue
            m=name2node[r.next_reference_name]
            if m==c: same+=1; continue
            oth+=1; L[m]+=1
            if r.mapping_quality>=20: L20[m]+=1
        b9=sum(v for k,v in L.items() if FINALBIN.get(k)=='bin.009'); b920=sum(v for k,v in L20.items() if FINALBIN.get(k)=='bin.009')
        hz=sum(v for k,v in L.items() if k in HZS)
        top='; '.join(f'{k}({v};{FINALBIN.get(k,"unbinned")})' for k,v in L.most_common(6))
        w.write('\t'.join(map(str,[c,INFO[c]['len'],n,un,same,oth,b9,b920,hz,top]))+'\n')
        for k,v in L.most_common(): wd.write(f'{c}\t{k}\t{FINALBIN.get(k,"unbinned")}\t{v}\t{L20[k]}\n')
