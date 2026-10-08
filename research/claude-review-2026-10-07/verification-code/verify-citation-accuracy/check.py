# Exact-integer check in the geometric representation of E_n (manuscript numbering:
# chain 0-1-...-(n-2), node n-1 attached to node 2).
from fractions import Fraction
def En_adj(n):
    adj={i:set() for i in range(n)}
    for i in range(n-2):
        adj[i].add(i+1); adj[i+1].add(i)
    adj[n-1].add(2); adj[2].add(n-1)
    return adj
def act(adj,n,word):
    # columns a_j = w(alpha_j) as integer vectors; right-multiply by s_i
    cols=[[1 if k==j else 0 for k in range(n)] for j in range(n)]
    for i in word:
        ai=cols[i]
        new=[c[:] for c in cols]
        new[i]=[-v for v in ai]
        for j in adj[i]:
            new[j]=[cols[j][k]+ai[k] for k in range(n)]
        cols=new
    return tuple(tuple(c) for c in cols)
def length(adj,n,w):
    # count positive roots sent negative via descent stripping
    cols=[list(c) for c in w]; L=0
    def neg(c): return any(v<0 for v in c)
    while True:
        d=[i for i in range(n) if neg(cols[i])]
        if not d: return L
        i=d[0]; L+=1
        ai=cols[i]; new=[c[:] for c in cols]; new[i]=[-v for v in ai]
        for j in adj[i]: new[j]=[cols[j][k]+ai[k] for k in range(n)]
        cols=new
def same(adj,n,w1,w2): return act(adj,n,w1)==act(adj,n,w2)
def digits(s): return [int(c) for c in s]

# Gern's D_m labels: s1,s2 attached to s3; chain s3-s4-...-sm.
# Gern's w_4 = s1 s2 s4 s3 s1 s2 s4 (Example 2.1.2).
gern_w4=[1,2,4,3,1,2,4]
def interval(i,j):
    # Gern Def 2.3.1 in labels 1..m
    if i>=2 and j>=i: return list(range(i,j+1))          # [i,j]=s_i...s_j
    if i==0: return [1,2]+list(range(3,j+1))             # [0,j]=s1 s2 s3..s_j
    if i==1: return [1]+list(range(3,j+1))               # [1,i]
    if j<i: return interval(j,i)[::-1]                   # [j,i]=[i,j]^-1 when j<i
    raise ValueError
def gern_wn(n):  # Lemma 2.3.4, n even
    k=n//2-2
    w=[]
    for j in range(2,n+1,2): w+=interval(j,0)
    # [n-k, n-2k] ... [n-1, n-2] [n,n]
    for t in range(k,-1,-1):
        w+=interval(n-t, n-2*t)
    return w

results=[]
# E_6: terminal 1325213, parabolic {1,2,3,5}: center 2 <-> s3; leaves 1,3,5 <-> s1,s2,s4
for n,word,phi in [(6,'1325213',{1:1,2:3,3:2,4:5}),
                   (7,'1326213',{1:1,2:3,3:2,4:6}),
                   (8,'1327213',{1:1,2:3,3:2,4:7})]:
    adj=En_adj(n); w=digits(word)
    img=[phi[g] for g in gern_w4]
    results.append((n,word,'w4 image',''.join(map(str,img)),same(adj,n,w,img),length(adj,n,act(adj,n,w))))
# D4 x A1 rows: w4 * commuting generator
for n,word,phi,u in [(7,'13256213',{1:1,2:3,3:2,4:6},5),
                     (8,'13257213',{1:1,2:3,3:2,4:7},5),
                     (8,'61327213',{1:1,2:3,3:2,4:7},6)]:
    adj=En_adj(n); w=digits(word)
    img=[phi[g] for g in gern_w4]+[u]
    results.append((n,word,'w4*u image',''.join(map(str,img)),same(adj,n,w,img),length(adj,n,act(adj,n,w))))
# D6 rows: Gern labels 1..6 -> 1, n-1, 2, 3, 4, 5  (manuscript line 388-389)
w6=gern_wn(6)
for n,word in [(7,'132543621324356'),(8,'132543721324357')]:
    adj=En_adj(n); phi={1:1,2:n-1,3:2,4:3,5:4,6:5}
    img=[phi[g] for g in w6]
    results.append((n,word,'w6 image',''.join(map(str,img)),same(adj,n,digits(word),img),length(adj,n,act(adj,n,digits(word)))))
for r in results: print(r)
print('len(gern w6 word)=',len(w6), w6)
