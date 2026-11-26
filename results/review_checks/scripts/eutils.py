import sys, urllib.request, urllib.parse, time, json
B="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
def get(ep, **p):
    p.setdefault("tool","review_checks"); 
    u=B+ep+"?"+urllib.parse.urlencode(p)
    for i in range(4):
        try:
            return urllib.request.urlopen(u, timeout=120).read().decode()
        except Exception as e:
            time.sleep(2+i*3)
    raise RuntimeError(u)
def efetch(ids, rettype="fasta", db="protein"):
    time.sleep(0.4); return get("efetch.fcgi", db=db, id=",".join(ids), rettype=rettype, retmode="text")
def esearch(term, db="protein", retmax=20):
    time.sleep(0.4); return json.loads(get("esearch.fcgi", db=db, term=term, retmax=retmax, retmode="json"))["esearchresult"]["idlist"]
def esummary(ids, db="protein"):
    time.sleep(0.4); r=json.loads(get("esummary.fcgi", db=db, id=",".join(ids), retmode="json"))["result"]
    return [r[u] for u in r["uids"]]
if __name__=="__main__":
    cmd=sys.argv[1]
    if cmd=="fetch": print(efetch(sys.argv[2].split(","), *(sys.argv[3:4] or ["fasta"])), end="")
    elif cmd=="search":
        db=sys.argv[3] if len(sys.argv)>3 else "protein"
        ids=esearch(sys.argv[2], db=db, retmax=int(sys.argv[4]) if len(sys.argv)>4 else 20)
        if ids:
            for d in esummary(ids, db=db): print(d.get("accessionversion",d.get("caption")), d.get("slen",""), d.get("title","")[:150], sep="\t")
