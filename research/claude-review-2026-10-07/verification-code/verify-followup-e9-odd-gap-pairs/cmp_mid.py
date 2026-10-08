"""Compare my vkl (DUMP mode) against the authors' independent targeted.py on random
words in E6, E9, D4, A4. For each target W, targeted.py gives P(x,W) for all x in [e,W];
vkl gives P(u w_J, W) for u in W^J. Since P(x,W)=P(x w_J... ) we compare against the
coset-maximal x = u*w_J in the authors' data, and also check the coset-constancy
P(x,W)=P(x t,W) for t in J in the authors' data.
"""
import sys, random, subprocess
sys.path.insert(0,'/Users/ic/workspace/kl/research/en_families')
import targeted
HERE='/private/tmp/claude-501/-Users-ic-workspace-kl/45cb4eb2-517e-46a6-a21d-77de64370e52/scratchpad/verify-followup-e9-odd-gap-pairs/'
def pstr(p):
    s=' + '.join((str(c) if i==0 else '%dq^%d'%(c,i)) for i,c in enumerate(p) if c)
    return s if s else '0'
def run(typ,g,n,word):
    W=g.elt(word)
    J=list(g.desc(W))
    # w_J
    wJ=g.e
    # longest element of W_J via BFS
    seen={wJ};q=[wJ]
    for x in q:
        for s in J:
            y=g.right(x,s)
            if y not in seen: seen.add(y);q.append(y)
    wJ=max(seen,key=g.length)
    out=subprocess.run([HERE+'vkl',typ]+[str(s) for s in word],capture_output=True,text=True,env={'DUMP':'1'}).stdout
    mine={}
    for line in out.split('\n'):
        if line.startswith('DUMP u='):
            uw,poly=line[7:].split(' : ')
            mine[tuple(int(c) for c in uw)]=poly
    theirs={}
    for x in g.lower(W):
        theirs[x]=pstr(g.kl(x,W))
    # check: every x in lower(W): P(x,W) equals mine at u = min rep of x's coset
    cnt=0;bad=0
    for x,p in theirs.items():
        # min rep: strip right descents in J
        u=x
        changed=True
        while changed:
            changed=False
            for t in J:
                if t in g.desc(u): u=g.right(u,t);changed=True
        uw=g.word(u)
        if uw not in mine:
            print('missing',uw);bad+=1;continue
        if mine[uw]!=p:
            print('MISMATCH',g.word(x),p,mine[uw]);bad+=1
        cnt+=1
    assert len(mine)*len(seen)==len(theirs),(len(mine),len(seen),len(theirs))
    return cnt,bad
random.seed(99)
total=0
for typ,n,g in [('A4',4,targeted.Coxeter(4,[(0,1),(1,2),(2,3)])),('D4',4,targeted.Coxeter(4,[(0,1),(1,2),(1,3)])),('E6',6,targeted.en(6)),('E9',9,targeted.en(9))]:
    for trial in range(3):
        L=random.randint(9,10 if typ=='E9' else 11)
        word=[];w=g.e
        while len(word)<L:
            s=random.randrange(n)
            if s in g.desc(w):continue
            w=g.right(w,s);word.append(s)
        cnt,bad=run(typ,g,n,word)
        total+=cnt
        print(typ,word,'pairs',cnt,'bad',bad,flush=True)
        if bad: sys.exit(1)
print('TOTAL pairs compared',total,'all equal')
