"""Exact all-parameter root certificates for terminal reflections in E10/E13.

Run: uv run python research/en_families/indefinite_positive_families.py
Only integral arithmetic and finite root reductions are used.  The KL
conclusion is a separate deduction from Green's maximum-independent FC lemma.
"""
from pathlib import Path
import json
import sys
from targeted import en

sys.setrecursionlimit(10000)

DATA = [
    dict(rank=10, a=1, beta=(3,7,10,9,7,6,4,3,1,6), I=(1,3,5,7,9),
         gamma=(1,2,3,2,2,1,1,0,0,1),
         matching=((0,1),(2,9),(3,4),(5,6),(7,8))),
    dict(rank=13, a=2, beta=(2,3,5,5,4,4,3,3,2,2,1,1,3),
         I=(0,3,5,7,9,11,12),
         gamma=(1,2,3,2,2,1,1,0,0,0,0,0,1),
         matching=((0,1),(2,12),(3,4),(5,6),(7,8),(9,10))),
]


def certificate(data):
    n=data['rank']; g=en(n); beta=data['beta']; gamma=data['gamma']; I=data['I']
    delta=tuple([2,4,6,5,4,3,2,1]+[0]*(n-9)+[3])

    def pairing(a,b):
        return sum(a[i]*(2*b[i]-sum(b[j] for j in g.adj[i])) for i in range(n))

    def mvec(a):
        return tuple(pairing(e,a) for e in g.e)

    def add(a,b):
        return tuple(x+y for x,y in zip(a,b))

    def sub(a,b):
        return tuple(x-y for x,y in zip(a,b))

    def reflect(root,a):
        m=pairing(root,a)
        return tuple(x-m*y for x,y in zip(a,root))

    def root_certificate(a):
        initial=a; steps=[]
        assert min(a)>=0 and pairing(a,a)==2
        while sum(a)>1:
            s=next(i for i in range(n) if pairing(g.e[i],a)>0)
            b=reflect(g.e[s],a)
            assert min(b)>=0 and sum(b)<sum(a)
            steps.append(s);a=b
        assert a in g.e
        simple=g.e.index(a)
        word=tuple(steps)+(simple,)+tuple(reversed(steps))
        matrix=tuple(reflect(initial,e) for e in g.e)
        assert g.elt(word)==matrix
        return dict(reduction=steps,simple_root=simple,reflection_word=word)

    def T(a):
        return reflect(gamma,reflect(add(gamma,delta),a))

    certs={name:root_certificate(root) for name,root in
           [('beta',beta),('gamma',gamma),('gamma_plus_delta',add(gamma,delta))]}
    tword=tuple(certs['gamma']['reflection_word'])+tuple(certs['gamma_plus_delta']['reflection_word'])
    tmatrix=tuple(T(e) for e in g.e)
    assert g.elt(tword)==tmatrix
    a=data['a']
    assert pairing(beta,gamma)==-3 and pairing(beta,delta)==-a
    for e in g.e:
        # This formula holds in the ambient indefinite root space, where
        # delta is not radical.
        assert T(e)==tuple(x+pairing(e,delta)*y-(pairing(e,gamma)+pairing(e,delta))*z
                            for x,y,z in zip(e,gamma,delta))
    v=sub(T(beta),beta);u=sub(T(v),v)
    assert v==tuple(-a*x+(a+3)*y for x,y in zip(gamma,delta))
    assert u==tuple(2*a*x for x in delta)
    # These three identities give T^k beta=beta+k*v+binom(k,2)*u by induction.
    assert T(beta)==add(beta,v) and T(v)==add(v,u) and T(u)==u
    assert pairing(delta,delta)==0 and pairing(gamma,delta)==0
    assert min(beta)>0 and min(v)>=0 and min(u)>=0 and any(v)
    coefficients=[mvec(a) for a in (beta,v,u)]
    for degree,m in enumerate(coefficients):
        assert all(m[i]>(0 if degree==0 else -1) for i in I)
        assert all(m[i]<=0 for i in range(n) if i not in I)
        assert all(m[s]+m[t]<=0 for s in I for t in g.adj[s])
    assert all(t not in I for s in I for t in g.adj[s])
    matching=data['matching'];vertices=[i for edge in matching for i in edge]
    assert len(vertices)==len(set(vertices)) and all(b in g.adj[a] for a,b in matching)
    assert len(I)==n-len(matching) and len(I)%2==1
    # Maximum size follows from this matching: an independent set uses at
    # most one endpoint of each matching edge plus every unmatched vertex.
    checks=[]
    for k in range(5):
        bk=tuple(x+k*y+k*(k-1)//2*z for x,y,z in zip(beta,v,u))
        mk=mvec(bk);w=tuple(reflect(bk,e) for e in g.e)
        assert pairing(bk,bk)==2 and g.inv(w)==w
        assert g.terminal(w) and g.desc(w)==I
        reduced=g.word(w)
        assert set(reduced)==set(range(n)) and len(reduced)%2==1
        assert g.leq(g.elt(I),w)
        checks.append(dict(k=k,beta=bk,pairings=mk,height=sum(bk),length=len(reduced),
                           descents=I,reduced_word=reduced,interval_rank=len(reduced)-len(I)))
        for fun in (g.right,g.desc,g.word,g.inv,g.leq):fun.cache_clear()
    return dict(**data,delta=delta,v=v,u=u,pairing_beta=coefficients[0],
                pairing_v=coefficients[1],pairing_u=coefficients[2],
                roots=certs,translation_word=g.word(tmatrix),
                translation_matrix_columns=tmatrix,
                polynomial_identity='T(beta)=beta+v; T(v)=v+u; T(u)=u',
                all_parameter_formula='beta_k=beta+k*v+binom(k,2)*u, k>=0',
                simplified_formula='beta_k=beta-a*k*gamma+(a*k*k+3*k)*delta, k>=0',
                checks=checks)


if __name__=='__main__':
    out=[certificate(data) for data in DATA]
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    for row in out:
        print(json.dumps({k:row[k] for k in ('rank','beta','I','delta','gamma','v','u',
                                          'pairing_beta','pairing_v','pairing_u')}))
        print('VERIFIED symbolic all-k identities, signs, maximum I, real roots; sample lengths',
              [x['length'] for x in row['checks']])
