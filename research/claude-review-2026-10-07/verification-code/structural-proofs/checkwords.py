import itertools
def setup(n):
    adj={i:set() for i in range(n)}
    for i in range(n-2): adj[i].add(i+1); adj[i+1].add(i)
    adj[n-1].add(2); adj[2].add(n-1)
    return adj
def rmul(w,i,adj,n):
    cols=list(w); ai=cols[i]; cols[i]=tuple(-c for c in ai)
    for j in adj[i]: cols[j]=tuple(a+b for a,b in zip(w[j],ai))
    return tuple(cols)
def elt(word,n,adj):
    w=tuple(tuple(1 if i==j else 0 for i in range(n)) for j in range(n))
    for s in word: w=rmul(w,s,adj,n)
    return w
def neg(c): return any(x<0 for x in c)
def R(w,n): return {j for j in range(n) if neg(w[j])}
def redlen(w,n,adj):
    k=0
    while R(w,n):
        w=rmul(w,min(R(w,n)),adj,n); k+=1
    return k
def inverse(w,n,adj):
    word=[]
    while R(w,n):
        s=min(R(w,n)); w=rmul(w,s,adj,n); word.append(s)
    return elt(word,n,adj)
def isFC(word,n,adj):
    # heap check: FC iff reduced word has no sts pattern reachable by commutations: use Stembridge criterion via all commutation-equivalent words is big; use the (left descent) BFS over right factors from checkwords
    pass
pairs={6:[("3125231","1325213")],7:[("3126231","1326213"),("31265231","13256213"),("653423126345231","132543621324356"),("6534231263452341230126345231","1325436210321432543621324356")]}
for n,lst in pairs.items():
    adj=setup(n)
    for mine,tab in lst:
        a=elt([int(c) for c in mine],n,adj); b=elt([int(c) for c in tab],n,adj)
        print(n,mine,tab,"same element:",a==b,"len",redlen(b,n,adj))
# E8 table rows: check lengths, descents, support, and that the x is a prefix (x<=b) and descents equal supp(x)
n=8; adj=setup(n)
rows=[("1327213","137"),("13257213","1357"),("61327213","1367"),("132543721324357","1357"),("1325437210321432543721324357","1357"),("7534231270123456210321432"+"5437210321432543721324357","1357")]
for bw,xw in rows:
    b=elt([int(c) for c in bw],n,adj); 
    L=R(inverse(b,n,adj),n); Rr=R(b,n)
    print("E8",bw[:20],"len",redlen(b,n,adj),"reduced word given:",redlen(b,n,adj)==len(bw),"L",sorted(L),"R",sorted(Rr),"supp x",sorted(int(c) for c in xw), "x subset supp(b):", set(int(c) for c in xw)<=set(int(c) for c in bw))
