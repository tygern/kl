"""An exact infinite affine-E8 family of full-support Gern-terminal elements.
Only 240 finite roots are enumerated. No affine group enumeration.
"""
from targeted import en
from pathlib import Path
import json,sys,math
sys.setrecursionlimit(10000)
g=en(9)
delta=(2,4,6,5,4,3,2,1,3)
beta=(1,2,3,3,2,2,1,1,2)
def pairing(a,b):return sum(a[i]*(2*b[i]-sum(b[j] for j in g.adj[i])) for i in range(9))
def add(a,b,k=1):return tuple(x+k*y for x,y in zip(a,b))
def action(w,a):return tuple(sum(w[j][i]*a[j] for j in range(9)) for i in range(9))
def root_reflection(a):
    return tuple(add(g.e[j],a,-pairing(g.e[j],a)) for j in range(9))
def reflect(a,s):
    out=list(a);out[s]-=pairing(g.e[s],a);return tuple(out)
def reflection_word(a):
    simple=[i for i in range(9) if a==g.e[i]]
    if simple:return (simple[0],)
    s=next(i for i in range(9) if pairing(g.e[i],a)>0)
    b=reflect(a,s);assert min(b)>=0 and sum(b)<sum(a)
    return (s,)+reflection_word(b)+(s,)

def finite_roots():
    inds=[i for i in range(9) if i!=7]
    roots={g.e[i] for i in inds};todo=list(roots)
    for a in todo:
        for s in inds:
            b=reflect(a,s)
            if b not in roots:roots.add(b);todo.append(b)
    return sorted(roots)

m=tuple(pairing(g.e[j],beta) for j in range(9));assert m==(0,0,-1,1,-1,1,-1,1,1)
assert all(pairing(delta,g.e[j])==0 for j in range(9)) and pairing(beta,beta)==2
# Root membership and an explicit conjugating translation certificate.
a=add(g.e[3],delta);rword=reflection_word(a);assert g.elt(rword)==root_reflection(a)
translation_word=rword+(3,);T=g.elt(translation_word)
assert action(T,beta)==add(beta,delta) and action(T,delta)==delta
bword=reflection_word(beta);assert g.elt(bword)==root_reflection(beta)
# Finite-root counting gives the affine inversion count symbolically.
b0=add(beta,delta,-1);roots=finite_roots();assert len(roots)==240
slope=offset=0;counts={}
for a in roots:
    val=pairing(a,b0)
    if val<=0:continue
    image=add(a,b0,-val);assert image in roots
    nmin=int(any(v<0 for v in a))
    image_negative=int(any(v<0 for v in image))
    # n>=nmin and n-val*(k+1)<0, plus the possible n=... boundary.
    slope+=val;offset+=-nmin+image_negative
    key=(val,nmin,image_negative);counts[key]=counts.get(key,0)+1
assert slope==58 and offset==-25
checks=[]
for k in range(11):
    bk=add(beta,delta,k);w=root_reflection(bk)
    assert g.length(w)==33+58*k and g.terminal(w)
    assert g.elt(g.word(w))==w
    row=dict(k=k,beta=bk,length=g.length(w),desc=g.desc(w),terminal=g.terminal(w),extensions=[])
    for s in (0,1):
        v=g.right(w,s)
        assert g.length(v)==34+58*k and g.terminal(v)
        x=g.elt(tuple(g.desc(v)));assert g.fc(x) and g.leq(x,v)
        row['extensions'].append(dict(generator=s,length=g.length(v),desc=g.desc(v),bottom_word=g.desc(v),bottom_length=g.length(x),interval_rank=g.length(v)-g.length(x)))
    checks.append(row)
    g.right.cache_clear();g.desc.cache_clear();g.word.cache_clear();g.inv.cache_clear()
out=dict(delta=delta,beta=beta,cartan_pairing=m,base_reflection_word=g.word(root_reflection(beta)),
         translation_word=translation_word,translation_reduced_word=g.word(T),
         finite_root_count=len(roots),length_slope=slope,length_offset_for_k_plus_1=offset,
         inversion_root_classes=[dict(pairing=key[0],minimum_affine_level=key[1],image_negative=key[2],count=value) for key,value in sorted(counts.items())],checks=checks)
Path(__file__).with_name('affine_reflection_family.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
