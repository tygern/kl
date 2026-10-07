"""Exact symbolic certificate for a full-support affine-E8 terminal family.
The polynomial conclusion uses the separately cited FC-cell uniqueness lemma;
no KL polynomial computation or whole affine group enumeration is performed.
"""
from targeted import en
from pathlib import Path
import json,sys
sys.setrecursionlimit(10000)
g=en(9);delta=(2,4,6,5,4,3,2,1,3);gamma=(1,1,1,1,1,0,0,0,0)
I=(1,3,5,7,8);d=(2,-1,1,-1,1,-1,1,-1,-1)
base_word=(3,1,2,8,7,5,6,4,5,2,3,4,1,2,3,0,1,2,8,5,6,7,4,5,2,3,1)
b=g.elt(base_word)
def pair(a,z):return sum(a[i]*(2*z[i]-sum(z[j] for j in g.adj[i])) for i in range(9))
def apply(w,a):return tuple(sum(w[j][i]*a[j] for j in range(9)) for i in range(9))
def mul(w,v):return tuple(apply(w,col) for col in v)
def refl(a):return tuple(tuple(int(i==j)-pair(g.e[j],a)*a[i] for i in range(9)) for j in range(9))
def reflect(a,s):
 z=list(a);z[s]-=pair(g.e[s],a);return tuple(z)
def rootword(a):
 if a in g.e:return (g.e.index(a),)
 s=next(i for i in range(9) if pair(g.e[i],a)>0);z=reflect(a,s)
 assert min(z)>=0 and sum(z)<sum(a)
 return (s,)+rootword(z)+(s,)
def finite_roots():
 roots=set(g.e[i] for i in range(9) if i!=7);todo=list(roots)
 for a in todo:
  for s in range(9):
   if s==7:continue
   z=reflect(a,s)
   if z not in roots:roots.add(z);todo.append(z)
 return sorted(roots)
assert g.length(b)==27 and g.inv(b)==b and g.terminal(b) and g.desc(b)==I
assert all(pair(delta,g.e[j])==0 for j in range(9)) and pair(gamma,gamma)==2
v=tuple(x-y for x,y in zip(gamma,apply(b,gamma)));assert tuple(pair(g.e[j],v) for j in range(9))==d
# Translation T acts on roots by a -> a - (a,gamma)delta.
gd=tuple(x+y for x,y in zip(gamma,delta))
tword=rootword(gamma)+rootword(gd);T=g.elt(tword)
Tformula=tuple(tuple(int(i==j)-pair(g.e[j],gamma)*delta[i] for i in range(9)) for j in range(9))
assert T==Tformula and g.length(T)==58
# Symbolic all-k terminal sign inequalities.
assert all(d[s]<=0 for s in I) and all(d[s]>=0 for s in range(9) if s not in I)
assert all(d[s]+d[t]>=0 for s in I for t in g.adj[s])
assert sum(delta[i]*d[i] for i in range(9))==0
assert tuple(sum(d[i]*b[j][i] for i in range(9)) for j in range(9))==tuple(-x for x in d)
# Exact all-k length via affine-root inversion counting (optional for parity).
roots=finite_roots();assert len(roots)==240
counts=[];intercept=slope=0
for a in roots:
 image=apply(b,a);level=image[7];eta=tuple(x-level*z for x,z in zip(image,delta));assert eta in roots
 nmin=int(any(x<0 for x in a));negative=int(any(x<0 for x in eta))
 c=-level-nmin+negative;sl=-sum(a[j]*d[j] for j in range(9))
 # Term is max(0,c+k*sl); the asserted signs prove it linear for all k>=0.
 if sl>0:assert c>=0;intercept+=c;slope+=sl
 elif sl<0:assert c<=0
 else:intercept+=max(c,0)
 counts.append(dict(root=a,image_finite=eta,image_level=level,minimum_level=nmin,linear_constant=c,linear_slope=sl))
assert (intercept,slope)==(27,92)
x=g.elt(I);checks=[]
for k in range(6):
 wk=tuple(tuple(b[j][i]+k*d[j]*delta[i] for i in range(9)) for j in range(9))
 assert g.length(wk)==27+92*k and g.terminal(wk) and g.desc(wk)==I and g.inv(wk)==wk
 assert g.leq(x,wk)
 if k==1:assert wk==mul(mul(T,b),g.inv(T))
 checks.append(dict(k=k,length=g.length(wk),descents=g.desc(wk),bottom_length=5,interval_rank=g.length(wk)-5,full_support=len(set(g.word(wk)))==9,reduced_word=g.word(wk)))
 g.right.cache_clear();g.desc.cache_clear();g.word.cache_clear();g.inv.cache_clear()
out=dict(numbering='Chain 0-1-2-3-4-5-6-7, node8 attached to2; affine node7.',delta=delta,gamma=gamma,I=I,
         base_word=base_word,base_matrix_columns=b,conjugating_translation_word=g.word(T),column_slopes=d,shift_vector=v,
         formula='b_k = T^k b_0 T^(-k); column_j(b_k)=column_j(b_0)+k*d_j*delta',
         exact_length_intercept=intercept,exact_length_slope=slope,
         finite_root_count=len(roots),inversion_count_certificate=counts,checks=checks)
Path(__file__).with_name('positive_affine_family.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ('inversion_count_certificate','checks')},indent=2))
print('VERIFIED',len(checks),'concrete elements and symbolic all-k terminal/length conditions')
