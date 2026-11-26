import sys; sys.path.insert(0,__file__.rsplit('/',1)[0]); from hits import *
H=load_tbn(); BP=load_bp()
GENES=[('hzsA','WP_099325812.1'),('hzsB','WP_164994616.1'),('hzsC','WP_478210066.1'),('hdh','Q1PW30.2'),
 ('nxrA (Nitrospira)','YEH58028.1'),('nxrA (anammox, kustd1700)','Q1PZD8'),('nxrB (Nitrospira)','YEH58030.1'),
 ('amoA','Q04507.2'),('amoB','Q04508.2'),('amoC (full-length)','WP_011111554.1'),('amoC (WP_234989093.1, partial ref)','WP_234989093.1'),
 ('amoA comammox clade A','WP_062484767.1'),('amoA comammox clade B','TKS59187.1'),
 ('nirS (Kuenenia)','WP_164995172.1'),('nirS (Paracoccus)','Q51700.1'),('norB','BAA32546.1'),('norC','Q51662.3'),
 ('napA','WP_104493079.1'),('narG','WP_230999932.1'),('nosZ clade I (Paracoccus)','WP_230999923.1'),('nosZ clade II (Desulfobacillus BBO22106.1)','BBO22106.1'),('mcrA','WP_505736924.1')]
R=[]
for g,ref in GENES: R+=rows(g,ref,H,BP)
write(RC+'/06_tableS2/TableS2_all_hits_e1e-10.tsv',R)
write(RC+'/06_tableS2/TableS2_hits_id40.tsv',[r for r in R if r[10]>=40])
print(len(R),'hits;',sum(r[10]>=40 for r in R),'with id>=40')
import collections
c=collections.Counter((r[0],r[-1]) for r in R if r[10]>=40)
for k,v in sorted(c.items()): print(k,v)
