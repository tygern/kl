from fractions import Fraction
from itertools import product
import sympy as sp
# E8 labeling of the manuscript: chain 0-1-2-3-4-5-6, node 7 attached to 2
n=8
edges={(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(2,7)}
def adj(i,j): return (i,j) in edges or (j,i) in edges
A=sp.zeros(n,n)
for i in range(n):
    for j in range(n):
        A[i,j]=2 if i==j else (-1 if adj(i,j) else 0)
def simple(i):
    M=sp.eye(n)
    for j in range(n):
        M[i,j]-=A[i,j]   # s_i(alpha_j)=alpha_j - A[i,j] alpha_i  (columns = images)
    return M
S=[simple(i) for i in range(n)]
word="7534231270123456210321432"+"5437210321432543721324357"
print(len(word))
W=sp.eye(n)
for c in word:
    W=W*S[int(c)]
print("w^2==I:", W*W==sp.eye(n))
x=sp.symbols('x')
print("charpoly:", sp.factor((W-x*sp.eye(n)).det()))
print("dim(-1)-eig:", n-(W+sp.eye(n)).rank())
print("dim(+1)-eig:", n-(W-sp.eye(n)).rank())
# length via inversion count: positive roots sent negative
# enumerate E8 positive roots by closure
roots=set()
fr=[tuple(1 if j==i else 0 for j in range(n)) for i in range(n)]
roots=set(fr); frontier=list(fr)
while frontier:
    new=[]
    for r in frontier:
        v=sp.Matrix(r)
        for i in range(n):
            w=S[i]*v
            t=tuple(int(c) for c in w)
            if all(c>=0 for c in t) and t not in roots:
                roots.add(t); new.append(t)
    frontier=new
print("num pos roots", len(roots))
inv=0
for r in roots:
    w=W*sp.Matrix(r)
    if all(int(c)<=0 for c in w): inv+=1
print("length", inv)
# -1 eigenspace basis in root coords
ns=(W+sp.eye(n)).nullspace()
print("(-1)-eigenspace:", [list(v.T) for v in ns])
# find orthogonal root pair whose product is W
pos=sorted(roots, key=lambda r: sum(r))
def refl(r):
    v=sp.Matrix(r)
    return sp.eye(n)-v*(v.T*A)
for r in pos:
    if (W*sp.Matrix(r)) == -sp.Matrix(r):
        print("root in -1 space:", r, "height", sum(r))
hits=[]
cand=[r for r in pos if (W*sp.Matrix(r))==-sp.Matrix(r)]
for a in cand:
    for b in cand:
        if a<b and refl(a)*refl(b)==W:
            hits.append((a,b,(sp.Matrix(a).T*A*sp.Matrix(b))[0]))
print("pairs", hits)
