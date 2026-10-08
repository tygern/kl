"""Enumerate all FC elements of affine E8 (E9 numbering: chain 0-..-7, node 8 on 2),
then list those x with I subset R(x), I subset L(x) and x <= W (via coset projection into
the dumped quotient ideal of W' = W w_I). Exact integer matrices only."""
import sys
n=9
edges=[(i,i+1) for i in range(7)]+[(2,8)]
A=[[0]*n for _ in range(n)]
for i in range(n):A[i][i]=2
for a,b in edges:A[a][b]=A[b][a]=-1
adj=[[j for j in range(n) if A[i][j]==-1] for i in range(n)]
E=tuple(tuple(int(i==j) for i in range(n)) for j in range(n))  # columns: E[j] = x(alpha_j)
def right(x,s):
    out=list(x);out[s]=tuple(-a for a in x[s])
    for t in adj[s]:out[t]=tuple(a+b for a,b in zip(x[s],x[t]))
    return tuple(out)
def rdesc(x):return [j for j in range(n) if all(a<=0 for a in x[j])]
# BFS of FC elements layer by layer (Stembridge: w FC iff every right-weak prefix has commuting descents)
layers=[{E}]
allfc={E:()}
while True:
    prev=layers[-1];cur={}
    for x in prev:
        ds=set(rdesc(x))
        for s in range(n):
            if s in ds:continue
            y=right(x,s)
            if y in cur:continue
            dy=rdesc(y)
            if any(A[a][b]==-1 for a in dy for b in dy):continue
            if all(right(y,t) in prev for t in dy):
                cur[y]=allfc[x]+(s,)
    if not cur:break
    layers.append(set(cur));allfc.update(cur)
    print('length',len(layers)-1,'FC count',len(cur),'cumulative',len(allfc),flush=True)
print('total FC elements',len(allfc),'max length',len(layers)-1)
# inverse via reversed word
def elt(word):
    x=E
    for s in word:x=right(x,s)
    return x
I=[3,5,7,8]
W=elt([8,7,5,6,3,4,5,2,3,4,1,2,3,0,1,2,8,2,3,4,1,2,3,0,1,2,8,5,6,7,4,5,3])
ideal=set()
for line in open(sys.argv[1]):
    vals=list(map(int,line.strip().split(',')))
    # dumped as M[i][j] row-major: M[i][j]=coeff of alpha_i in x(alpha_j); our x[j][i]
    cols=tuple(tuple(vals[i*n+j] for i in range(n)) for j in range(n))
    ideal.add(cols)
print('ideal loaded',len(ideal))
def proj(x):  # strip right descents in I
    changed=True
    while changed:
        changed=False
        for t in I:
            if all(a<=0 for a in x[t]): x=right(x,t);changed=True
    return x
elig=[]
for x,word in allfc.items():
    if not set(I)<=set(rdesc(x)):continue
    xi=elt(reversed(word))
    if not set(I)<=set(rdesc(xi)):continue
    if proj(x) in ideal: elig.append((len(word),word))
elig.sort()
print('eligible FC bottoms for W (I subset R and L, x<=W):',len(elig))
for l,w in elig:print('  length',l,'word',w)
