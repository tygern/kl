import sys, random, subprocess, json
sys.path.insert(0,'/Users/ic/workspace/kl/research/en_families')
import targeted
def run(typ, g, n, word):
    w=g.elt(word)
    mine={}
    out=subprocess.run(['./klq','dump',typ]+[str(s) for s in word],capture_output=True,text=True).stdout
    for line in out.strip().split('\n'):
        key,poly=line.split(' ',1); mine[key]=poly
    theirs={}
    for x in g.lower(w):
        p=g.kl(x,w)
        cols=[]
        for j in range(n):
            c=list(x[j])+[0]*(n-len(x[j]))
            cols.append(','.join(str(v) for v in c)+';')
        key=''.join(cols)
        s=' + '.join((str(c) if i==0 else '%dq^%d'%(c,i)) for i,c in enumerate(p) if c)
        theirs[key]=s
    assert set(mine)==set(theirs), (len(mine),len(theirs))
    bad=[k for k in mine if mine[k]!=theirs[k]]
    return len(mine), bad
random.seed(5)
total=0
for typ,n,g in [('E6',6,targeted.en(6)),('E9',9,targeted.en(9)),('D4',4,targeted.Coxeter(4,[(0,1),(1,2),(1,3)]))]:
    for trial in range(6):
        L=random.randint(6,11 if typ=='E9' else 12)
        word=[]
        w=g.e
        while len(word)<L:
            s=random.randrange(n)
            if s in g.desc(w): continue
            w=g.right(w,s); word.append(s)
        cnt,bad=run(typ,g,n,word)
        total+=cnt
        print(typ,word,'pairs',cnt,'mismatches',len(bad))
        if bad: print('  BAD',bad[:3]); sys.exit(1)
print('total pairs compared',total)
