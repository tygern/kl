"""Save every KL dependency for the affine D4 star example, and bounded family checks."""
from targeted import star
from itertools import product
from math import comb
from collections import Counter
from pathlib import Path
import json

OUT=Path(__file__).resolve().parent

def certify(m,certificate=False,all_orders=False):
    g=star(m);A=tuple(range(1,m+1));x=g.elt(A);w=g.elt(A+(0,)+A)
    p=g.kl(x,w);expected=tuple(comb(m,k)-(comb(m,k-1) if k else 0) for k in range(m//2+1));assert p==expected
    lower=g.lower(w);interval=[z for z in lower if g.leq(x,z)]
    forms={}
    for statuses in product(range(4),repeat=m):
        L=tuple(i+1 for i,s in enumerate(statuses) if s&1)
        R=tuple(i+1 for i,s in enumerate(statuses) if s&2)
        z=g.elt(L+(0,)+R)
        assert z not in forms;forms[z]=(frozenset(L),frozenset(R))
        assert g.length(z)==len(L)+len(R)+1
        assert all((j in g.desc(z))==(j in R) for j in A)
        if z in interval:assert set(L)|set(R)==set(A)
    assert len(lower)==4**m+2**m and len(interval)==3**m+1
    assert set(g.eligible(w))==({x} if m>=2 else set())
    assert g.terminal(w)==(m>=3)
    if all_orders:
        for z,(L,R) in forms.items():
            for z2,(L2,R2) in forms.items():
                assert g.leq(z,z2)==(L<=L2 and R<=R2)
    for s in g.desc(w):assert g.with_descent(x,w,s)==p
    ranks=Counter(g.length(z)-m for z in interval)
    out=dict(leaves=m,top_word=A+(0,)+A,bottom_word=A,length_top=2*m+1,length_bottom=m,
             lower_size=len(lower),interval_size=len(interval),rank_vector=[ranks[i] for i in range(m+2)],polynomial=p,
             mu=p[m//2] if m%2==0 else 0,terminal=g.terminal(w),eligible_bottoms=len(g.eligible(w)),all_pair_orders_checked=all_orders)
    if certificate:
        # Materialize all correction records before fixing the element-ID dictionary.
        for a,b in list(g.kl_values):
            if a!=b and g.leq(a,b):g.corrections(b,g.desc(b)[0])
        elements=sorted(set(lower)|{z for pair in g.kl_values for z in pair})
        indices={z:i for i,z in enumerate(elements)}
        records=[]
        for (a,b),poly in g.kl_values.items():
            row=dict(x=indices[a],w=indices[b],polynomial=poly)
            if a!=b and g.leq(a,b):
                s=g.desc(b)[0];row['right_descent']=s
                row['correction_terms']=[(indices[z],d,mu) for z,d,mu in g.corrections(b,s) if g.leq(a,z)]
            records.append(row)
        data=dict(format_version=1,group=dict(rank=m+1,edges=g.edges),
             convention='Generator 0 is the center; generators 1,...,m are commuting leaves. Integral matrices are stored by columns. Polynomial coefficients ascend in q.',
             elements=[dict(matrix=z,word=g.word(z),length=g.length(z),R=g.desc(z)) for z in elements],
             root=[indices[x],indices[w]],records=records,
             top_recurrences=[dict(right_descent=s,correction_terms=[(indices[z],d,mu) for z,d,mu in g.corrections(w,s) if g.leq(x,z)]) for s in g.desc(w)],
             summary=out)
        (OUT/'affine_d4_kl_certificate.json').write_text(json.dumps(data,separators=(',',':'))+'\n')
    return out

if __name__=='__main__':
    out=[certify(m,certificate=m==4,all_orders=m<=4) for m in range(1,9)]
    (OUT/'star_family_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
