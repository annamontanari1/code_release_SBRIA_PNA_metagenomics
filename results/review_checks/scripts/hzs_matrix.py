import collections,os
os.chdir(os.path.expanduser('~/sbria_anammox_shotgun/results/review_checks/06_tableS2'))
# B. fulgida: rótulos do GenBank trocados; KKO20886.1 ('hzsB') = gamma (85,8% Q1Q0T3) -> hzsC; KKO20885.1 ('hzsC') = beta (80,4% Q1Q0T4) -> hzsB
# (nomenclatura Swiss-Prot/PDB 5C2V; 14_ncycdb/followup_21out, 21/10/2026)
sp={'GAN34140.1':('hzsA','B. sinica'),'GAN34138.1':('hzsB','B. sinica'),'GAN34139.1':('hzsC','B. sinica'),'WP_052564202.1':('hdh','B. sinica'),
'WP_080325085.1':('hzsA','B. sapporoensis'),'WP_080325083.1':('hzsB','B. sapporoensis'),'WP_080325084.1':('hzsC','B. sapporoensis'),'OQD44929.1':('hdh','B. sapporoensis'),
'KKO20072.1':('hzsA','B. fulgida'),'KKO20886.1':('hzsC','B. fulgida'),'KKO20885.1':('hzsB','B. fulgida'),'KKO18440.1':('hdh','B. fulgida'),
'WP_099325812.1':('hzsA','Ca. Kuenenia'),'WP_164994616.1':('hzsB','Ca. Kuenenia'),'WP_478210066.1':('hzsC','Ca. Kuenenia'),'Q1PW30.2':('hdh','Ca. Kuenenia')}
SPP=['B. sinica','B. sapporoensis','B. fulgida','Ca. Kuenenia']
d=collections.defaultdict(dict); qlen={}
for l in open('hzs_orfs_vs_anammox_panel_long.tsv'):
    p=l.rstrip().split('\t'); g,s=sp[p[1]]; qlen[p[0]]=p[10]
    cov=100*(int(p[5])-int(p[4])+1)/int(p[10]); k=(p[0],g)
    if s not in d[k] or float(p[9])>d[k][s][2]: d[k][s]=(float(p[2]),round(cov,1),float(p[9]),p[1])
with open('hzs_orfs_vs_Brocadia_Kuenenia_matrix.tsv','w') as w:
    w.write('ORF\tORF_len_aa\tgene\t'+'\t'.join(f'{s} id% (ORFcov%)' for s in SPP)+'\tbest_match\n')
    for (o,g),v in sorted(d.items()):
        if max(x[0] for x in v.values())<40: continue
        best=max(v.items(),key=lambda x:x[1][2])[0]
        w.write(f'{o}\t{qlen[o]}\t{g}\t'+'\t'.join(f'{v[s][0]:.1f} ({v[s][1]})' if s in v else '-' for s in SPP)+f'\t{best}\n')
