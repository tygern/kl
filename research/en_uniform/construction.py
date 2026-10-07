"""Counterchecks for the rank-uniform E_(4r+1) construction, r>=3.

The all-r proof is in construction.txt; finite samples here do not replace it.
No group-ball, Bruhat-interval, or FC-catalogue enumeration is used.
"""
from pathlib import Path
import argparse
import json
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'en_families'))
from targeted import en


def seed(r):
    return (r-1,r,2*r-1)+tuple((4*r+1-j)//2 for j in range(3,4*r))+(r,)


def verify(r,ks):
    n=4*r+1;N=4*r;g=en(n);b=seed(r);a=2*r-4
    I=(0,)+tuple(range(3,4*r,2))+(4*r,)
    delta=(2,4,6,5,4,3,2,1)+(0,)*(n-9)+(3,)
    gamma=(1,2,3,2,2,1,1,0)+(0,)*(n-9)+(1,)

    def pair(x,y):
        return sum(x[i]*(2*y[i]-sum(y[j] for j in g.adj[i])) for i in range(n))

    def reflect(root,x):
        m=pair(root,x)
        return tuple(z-m*t for z,t in zip(x,root))

    def dcoords_twice(c):
        # Twice the finite-D component in alpha0=h+lambda coordinates.
        return ((-c[0]+2*c[1]-2*c[N]),
                (c[0]-2*c[1]+2*c[2]-2*c[N]))+tuple(
                    c[0]-2*c[j-1]+2*c[j] for j in range(3,N))+(c[0]-2*c[N-1],)

    after_s0=reflect(g.e[0],b)
    assert after_s0[0]==1 and after_s0[1:]==b[1:]
    lambdatwice=(-1,)+(1,)*(N-1)
    signs=tuple(-1 if j%2==0 else 1 for j in range(1,N+1))
    assert sum(x<0 for x in signs)==2*r  # even signed permutation, in W(D_N)
    assert dcoords_twice(after_s0)==tuple(x*y for x,y in zip(signs,lambdatwice))
    assert pair(b,b)==2 and min(b)>0
    assert pair(gamma,gamma)==2 and pair(delta,delta)==pair(gamma,delta)==0
    assert pair(b,gamma)==-r and pair(b,delta)==-a
    gd=tuple(x+y for x,y in zip(gamma,delta))
    gamma_steps=(2,1,0,4,3,2,1,6,5,4,3,2)
    gd_steps=gamma_steps+(N,2,1,0,3,2,1,4,3,2,5,4,3,6,5,4,7,6,5,N,2,1,0,3,2,1,4,3,2)
    for root,steps in ((gamma,gamma_steps),(gd,gd_steps)):
        for s in steps:
            old_height=sum(root);root=reflect(g.e[s],root)
            assert min(root)>=0 and sum(root)<old_height
        assert root==g.e[N]

    def T(x):return reflect(gamma,reflect(gd,x))

    for e in g.e:
        assert T(e)==tuple(x+pair(e,delta)*y-(pair(e,gamma)+pair(e,delta))*z
                           for x,y,z in zip(e,gamma,delta))
    v=tuple((a+r)*d-a*c for c,d in zip(gamma,delta))
    u=tuple(2*a*d for d in delta)
    assert min(v)>=0 and min(u)>=0 and any(v)
    assert T(b)==tuple(x+y for x,y in zip(b,v))
    assert T(v)==tuple(x+y for x,y in zip(v,u)) and T(u)==u
    coeffs=[tuple(pair(e,z) for e in g.e) for z in (b,v,u)]
    for degree,m in enumerate(coeffs):
        assert all(m[i]>=(1 if degree==0 else 0) for i in I)
        assert all(m[i]<=0 for i in range(n) if i not in I)
        assert all(m[s]+m[t]<=0 for s in I for t in g.adj[s])
    matching=((0,1),(2,N))+tuple((j,j+1) for j in range(3,N-1,2))
    assert len(matching)==2*r and len(set(i for e in matching for i in e))==4*r
    assert len(I)==2*r+1 and all(t not in I for s in I for t in g.adj[s])
    checks=[]
    for k in ks:
        root=tuple(z-a*k*c+(a*k*k+r*k)*d for z,c,d in zip(b,gamma,delta))
        expected=[r-2,2-r,-1-a*k]+[(-1)**(j+1) for j in range(3,N)]+[1+a*k]
        for j in range(3,8):expected[j]=(-1)**(j+1)*(1+a*k)
        expected[8]=-1-a*k*k-r*k
        assert min(root)>0 and pair(root,root)==2
        assert tuple(pair(e,root) for e in g.e)==tuple(expected)
        w=tuple(reflect(root,e) for e in g.e)
        assert g.desc(w)==I
        # Right terminality and involution imply left terminality without
        # reconstructing increasingly long reduced words.
        assert all(t not in g.desc(g.right(w,s)) for s in I for t in g.adj[s])
        assert all(reflect(root,reflect(root,e))==e for e in g.e)
        checks.append(dict(k=k,height=sum(root),largest_coordinate=max(root),
                           norm=2,descents=I,terminal=True))
    return dict(r=r,rank=n,seed=b,I=I,delta=delta,gamma=gamma,a=a,
                seed_s0_D_coordinates_twice=dcoords_twice(after_s0),
                even_sign_change_coordinates=list(range(2,N+1,2)),
                pairing_coefficients=coeffs,checks=checks)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--max-r',type=int,default=30)
    args=p.parse_args();assert args.max_r>=3
    rows=[verify(r,(0,1,2,10)) for r in range(3,args.max_r+1)]
    out=dict(scope='All-r proof is construction.txt; these are finite counterchecks',
             theorem_r_minimum=3,checked_r_maximum=args.max_r,
             parameter_formula='beta_(r,k)=beta_r-(2r-4)k*gamma+((2r-4)k^2+r*k)*delta',
             rows=rows)
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print(f'Verified r=3,...,{args.max_r}; ranks13,...,{4*args.max_r+1}; k=0,1,2,10.')
    print('Checked finite-D signed-coordinate seed certificate and symbolic translation identities.')
