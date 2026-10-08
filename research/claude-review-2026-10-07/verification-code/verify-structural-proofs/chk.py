# E8 in the paper's labeling: chain 0-1-2-3-4-5-6, node 7 attached to node 2.
n=8
edges={(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(2,7)}
def adj(i,j): return (i,j) in edges or (j,i) in edges
C=[[2 if i==j else (-1 if adj(i,j) else 0) for j in range(n)] for i in range(n)]
# vectors in root coordinates; s_i(v) = v - <v,alpha_i^vee> alpha_i, <v,alpha_i^vee> = sum_j v_j C[j][i]
def s(i,v):
    p=sum(v[j]*C[j][i] for j in range(n))
    w=list(v); w[i]-=p; return w
def act(word,v):
    for i in reversed(word): v=s(i,v)
    return v
word="7534231270123456210321432"+"5437210321432543721324357"
w=[int(c) for c in word]
# reducedness: w = s_{w0} s_{w1} ... ; word reduced iff each prefix s_{w0}..s_{w_{t-1}} applied... use: ell(w)=#positive roots sent negative.
# Enumerate positive roots of E8 by closure under simple reflections from simple roots.
simple=[[1 if j==i else 0 for j in range(n)] for i in range(n)]
roots=set(tuple(a) for a in simple); frontier=list(roots)
while frontier:
    new=[]
    for r in frontier:
        for i in range(n):
            q=tuple(s(i,list(r)))
            if q not in roots: roots.add(q); new.append(q)
    frontier=new
pos=[r for r in roots if all(c>=0 for c in r)]
print("number of roots",len(roots),"positive",len(pos))
def length(word):
    return sum(1 for r in pos if any(c<0 for c in act(word,list(r))))
print("word length",len(w),"Coxeter length",length(w),"reduced?",length(w)==len(w))
# matrix of w
M=[act(w,list(simple[j])) for j in range(n)]  # columns = images of alpha_j
M=[[M[j][i] for j in range(n)] for i in range(n)]
M2=[[sum(M[i][k]*M[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
print("involution?",M2==[[1 if i==j else 0 for j in range(n)] for i in range(n)])
# rank of M+I over Q via fraction-free elimination -> dim(-1)-eigenspace = n - rank(M+I)
from fractions import Fraction
A=[[Fraction(M[i][j]+(1 if i==j else 0)) for j in range(n)] for i in range(n)]
def rank(A):
    A=[row[:] for row in A]; r=0; m=len(A); k=len(A[0])
    for c in range(k):
        piv=next((i for i in range(r,m) if A[i][c]!=0),None)
        if piv is None: continue
        A[r],A[piv]=A[piv],A[r]
        for i in range(m):
            if i!=r and A[i][c]!=0:
                f=A[i][c]/A[r][c]; A[i]=[A[i][j]-f*A[r][j] for j in range(k)]
        r+=1
    return r
print("dim (-1)-eigenspace =",n-rank(A))
tr=sum(M[i][i] for i in range(n)); print("trace",tr,"(reflection would have trace",n-2,")")
# descent sets
L=[i for i in range(n) if any(c<0 for c in act(w+[i],list(simple[i])) ) ]  # placeholder
R=[i for i in range(n) if any(c<0 for c in act(w,list(simple[i])))]
Lset=[i for i in range(n) if any(c<0 for c in act([i]+w,[0]*n)) ]
# proper left descents: i in L(w) iff ell(s_i w)<ell(w)
Lset=[i for i in range(n) if length([i]+w)<len(w)]
print("right descents",R,"left descents",Lset)
# reflection lengths in finite Weyl group are odd: 2*ht-1
print("all reflection lengths odd; 50 is even ->", 50%2==0)
