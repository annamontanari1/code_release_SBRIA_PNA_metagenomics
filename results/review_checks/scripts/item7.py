import sys; sys.path.insert(0,__file__.rsplit('/',1)[0]); from base import *
import numpy as np
os.makedirs(RC+'/07_depth',exist_ok=True)
def wq(vals,w,qs):
    o=np.argsort(vals); v=np.array(vals)[o]; ww=np.array(w,float)[o]; c=np.cumsum(ww)-0.5*ww; c/=ww.sum()
    return [float(np.interp(q,c,v)) for q in qs]
def wpct(vals,w,x):
    v=np.array(vals); ww=np.array(w,float); return 100*ww[v<x].sum()/ww.sum()+50*ww[v==x].sum()/ww.sum()
BINS=['bin.004','bin.007','bin.009','bin.011','bin.016','bin.008','bin.017']
stats={}
with open(RC+'/07_depth/bin_depth_distribution_lenweighted.tsv','w') as w:
    w.write('bin\tn_contigs\ttotal_bp\tdepth_wmean\tdepth_wmedian\tdepth_wQ1\tdepth_wQ3\tdepth_wIQR\tdepth_wP5\tdepth_wP95\tdepth_min\tdepth_max\tGC_wmean\tGC_wsd\tGC_wP5\tGC_wP95\n')
    for b in BINS+[f'bin.{i:03d}' for i in range(1,26) if f'bin.{i:03d}' not in BINS]:
        cs=bin_contigs(b); d=[INFO[c]['depth'] for c in cs]; L=[INFO[c]['len'] for c in cs]; g=[INFO[c]['gc'] for c in cs]
        md,q1,q3,p5,p95=wq(d,L,[.5,.25,.75,.05,.95]); gm=np.average(g,weights=L); gsd=np.sqrt(np.average((np.array(g)-gm)**2,weights=L)); g5,g95=wq(g,L,[.05,.95])
        stats[b]=(d,L,g)
        w.write('\t'.join(map(str,[b,len(cs),sum(L),round(np.average(d,weights=L),2),round(md,2),round(q1,2),round(q3,2),round(q3-q1,2),round(p5,2),round(p95,2),round(min(d),2),round(max(d),2),round(gm,2),round(gsd,2),round(g5,2),round(g95,2)]))+'\n')
CAND=[('NODE_2432','nxrA (Nitrospira-type; added to bin.004_refined)','bin.004'),('NODE_719','nxrA fragment','bin.004'),
 ('NODE_2578','amoC (comammox-type)','bin.007'),('NODE_13511','amoA/amoB (comammox clade A)','bin.007'),('NODE_4141','nxrA','bin.007'),('NODE_348','nxr (log §11)','bin.007'),
 ('NODE_696','nirS','bin.009'),('NODE_4298','nxrA paralog','bin.009'),('NODE_1742','nxrA paralog','bin.009'),('NODE_11','hdh-like (HAO family)','bin.009'),
 ('NODE_5668','hzsA','bin.009'),('NODE_5312','hzsABC','bin.009'),('NODE_7170','hzsABC','bin.009'),('NODE_17695','hzsBC','bin.009'),('NODE_29475','hzsA','bin.009'),
 ('NODE_23806','hdh','bin.009'),('NODE_20914','hdh','bin.009'),('NODE_24842','hdh','bin.009'),('NODE_47442','hdh','bin.009'),
 ('NODE_46','narG; nosZ (bin.011)','bin.011'),('NODE_174','nirS','bin.011'),('NODE_5980','napA','bin.011'),
 ('NODE_929','napA/norB/norC/nirBD','bin.016'),('NODE_276','nirS','bin.016'),
 ('NODE_522','nirS','bin.008'),('NODE_2924','narG','bin.008'),('NODE_395','napA','bin.008'),
 ('NODE_6','napA/nirS','bin.017'),('NODE_57','norB/norC/nosZ','bin.017'),
 ('NODE_4969','amoCAB (Nitrosomonas)','none'),('NODE_43702','mcrA (Methanothrix)','none'),('NODE_5830','norB 42% (unbinned)','bin.008'),('NODE_25906','norC 59% (unbinned)','bin.008')]
with open(RC+'/07_depth/candidate_contig_positions.tsv','w') as w:
    w.write('contig\tgene(s)\tcompared_bin\tcontig_final_bin\tcontig_len\tdepth\tGC\tdepth_percentile_in_bin(len-weighted)\tGC_percentile_in_bin\tdepth/bin_wmedian\tinside_P5-P95_depth\tinside_P5-P95_GC\tbinners\n')
    for c,gn,b in CAND:
        I=INFO[c]
        if b=='none':
            w.write('\t'.join(map(str,[c,gn,'-',FINALBIN.get(c,'unbinned'),I['len'],round(I['depth'],2),I['gc'],'','','','','',binners(c)]))+'\n'); continue
        d,L,g=stats[b]; md,p5,p95=wq(d,L,[.5,.05,.95]); g5,g95=wq(g,L,[.05,.95])
        w.write('\t'.join(map(str,[c,gn,b,FINALBIN.get(c,'unbinned'),I['len'],round(I['depth'],2),I['gc'],round(wpct(d,L,I['depth']),1),round(wpct(g,L,I['gc']),1),round(I['depth']/md,2),p5<=I['depth']<=p95,g5<=I['gc']<=g95,binners(c)]))+'\n')
