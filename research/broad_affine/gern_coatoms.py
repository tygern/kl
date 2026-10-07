"""Count coatoms of Gern bad intervals by all signed transpositions."""
import sys,json
sys.path.insert(0,'computations')
from sparse_kl import SparseCoxeter
rows=[]
for n in range(4,25,2):
 g=SparseCoxeter('D',n)
 w=tuple((-1)**(n//2) if i==1 else i if i%2 else -(n+2-i) for i in range(1,n+1))
 x=tuple([-1,-2]+[i+1 if i%2 else i-1 for i in range(3,n+1)])
 coatoms=[];allcovers=0
 for i in range(n):
  for j in range(i+1,n):
   for sign in [1,-1]:
    z=list(w);z[i],z[j]=sign*w[j],sign*w[i];z=tuple(z)
    if g.length(z)==g.length(w)-1:
     allcovers+=1
     if g.leq(x,z):coatoms.append((i+1,j+1,sign))
 rows.append(dict(n=n,length=g.length(w),rank=g.length(w)-g.length(x),coatoms=len(coatoms),all_lower_covers=allcovers,cover_reflections=coatoms))
 print(rows[-1],flush=True)
with open('research/broad_affine/gern_coatoms.json','w') as f:json.dump(rows,f,indent=2)
