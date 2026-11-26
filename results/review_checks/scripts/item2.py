import sys; sys.path.insert(0,__file__.rsplit('/',1)[0]); from hits import *
os.makedirs(RC+'/02_amo_reflen',exist_ok=True)
H=load_tbn(); BP=load_bp()
w=open(RC+'/02_amo_reflen/amoB_amoC_coverage_fullrefs.tsv','w')
w.write('ref\tref_len_official\tref_len_in_amo_refs.fasta\tcontig\tcontig_len\tdepth\tGC\tfinal_bin\tbinners\tpident_bestHSP\taln_len\tqstart\tqend\tevalue\tbits\tref_cov_bestHSP_pct\tref_cov_union_pct\tblastp_ORF\tORF_len_aa\tORF_pid\tORF_qcov\tORF_scov\tstatus_cov70_id40\n')
old={'Q04508.2':70,'WP_234989093.1':70,'WP_011111554.1':'n.a.'}
for ref in ['Q04508.2','WP_234989093.1','WP_011111554.1']:
    for h in per_contig(H,ref):
        c=h['contig']; I=INFO[c]
        bp=[b for b in BP if b['ref']==ref and b['contig']==c]; b=max(bp,key=lambda x:x['bits']) if bp else None
        cov=h['qcov_union']; st='supported' if (h['pid']>=40 and cov>=70) else ('partial' if h['pid']>=40 else 'not supported (id<40)')
        w.write('\t'.join(map(str,[ref,REFLEN[ref],old[ref],c,I['len'],round(I['depth'],2),I['gc'],FINALBIN.get(c,'unbinned'),binners(c),h['pid'],h['aln'],h['qs'],h['qe'],h['ev'],h['bits'],round(h['qcov_best'],1),round(cov,1),
          b['orf'] if b else '',b['slen'] if b else '',b['pid'] if b else '',b['qcov'] if b else '',b['scov'] if b else '',st]))+'\n')
w.close()
