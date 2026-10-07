"""Independent exact audit of the rank-uniform E_(4r+1) family.

This imports neither the construction nor the family/FC engines.  The
all-rank proof is in referee.txt.  Symbolic assertions below check its
local identities in Z[r,k,q]; finite matrix checks are separate diagnostics.
"""
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


class P:
    """Small exact polynomial ring Q[r,k,q], for the printed identities."""
    def __init__(self, x=0):
        self.d = dict(x.d) if isinstance(x, P) else (
            {a: Fraction(b) for a, b in x.items() if b} if isinstance(x, dict)
            else ({(0, 0, 0): Fraction(x)} if x else {}))

    def __add__(self, x):
        d = dict(self.d)
        for a, b in P(x).d.items():
            d[a] = d.get(a, 0) + b
        return P(d)

    __radd__ = __add__

    def __neg__(self):
        return P({a: -b for a, b in self.d.items()})

    def __sub__(self, x):
        return self + -P(x)

    def __rsub__(self, x):
        return P(x) + -self

    def __mul__(self, x):
        d = {}
        for a, b in self.d.items():
            for c, e in P(x).d.items():
                f = tuple(v + w for v, w in zip(a, c))
                d[f] = d.get(f, 0) + b * e
        return P(d)

    __rmul__ = __mul__

    def __eq__(self, x):
        return self.d == P(x).d


def symbolic():
    r, k, q = [P({tuple(int(j == i) for j in range(3)): 1}) for i in range(3)]
    a = 2*r - 4
    # Pair cancellation on the long arm leaves exactly this seed norm.
    assert (r - 1)*(r - 2) + r*(2 - r) + r == 2
    # In the independently chosen D-coordinate realization, s_0 beta has
    # alpha_0 coefficient one and alternating finite coordinates -1/2,+1/2.
    assert r - 1 - (r - 2) == 1
    assert -(r - 1)*Fraction(1, 2) + (r - 2)*Fraction(1, 2) == Fraction(-1, 2)
    assert (3 - r)*Fraction(1, 2) + (r - 2)*Fraction(1, 2) == Fraction(1, 2)
    # Generic interior tail pairing, at even 2q and odd 2q+1.
    even, odd = 2*r - q, 2*r - q
    assert 2*even - (2*r - q + 1) - odd == -1
    assert 2*odd - even - (2*r - q - 1) == 1
    # The only nonconstant local section is supported in nodes 0..8,branch.
    delta = (2, 4, 6, 5, 4, 3, 2, 1, 0, 0)
    gamma = (1, 2, 3, 2, 2, 1, 1, 0, 0, 0)
    seed = [r-1, r] + [2*r-j//2 for j in range(2, 10)]
    beta = [seed[j] - a*k*gamma[j] + (a*k*k+r*k)*delta[j] for j in range(10)]
    branch = r-a*k+(a*k*k+r*k)*3
    pairings = [2*beta[0]-beta[1], 2*beta[1]-beta[0]-beta[2],
                2*beta[2]-beta[1]-beta[3]-branch]
    pairings += [2*beta[j]-beta[j-1]-beta[j+1] for j in range(3, 9)]
    expected = [r-2, 2-r, -1-a*k] + [(-1)**(j+1)*(1+a*k) for j in range(3, 8)]
    expected += [-1-a*k*k-r*k]
    assert pairings == expected
    assert 2*branch-beta[2] == 1+a*k
    m0 = [r-2, 2-r, -1, 1, -1, 1, -1, 1]
    assert sum(m0[j]*gamma[j] for j in range(8)) + 1 == -r
    assert sum(m0[j]*delta[j] for j in range(8)) + 3 == -a
    # Exceptional edge inequalities have nonpositive binomial coefficients.
    assert expected[7]+expected[8] == -r*k-a*k*(k-1)
    assert expected[8]+1 == -(a+r)*k-a*k*(k-1)
    assert beta[0] == r-1+4*k+(4*r-8)*k*k
    return {"coefficient_ring": "Q[r,k,q]", "all_assertions_pass": True,
            "edge_7_8_sum": "-r*k-(2*r-4)*k*(k-1)",
            "edge_8_9_sum": "-(3*r-4)*k-(2*r-4)*k*(k-1)"}


def check_rank(r):
    n, a = 4*r+1, 2*r-4
    edges = [(j, j+1) for j in range(n-2)] + [(2, n-1)]
    adj = [[] for _ in range(n)]
    for s, t in edges:
        adj[s].append(t); adj[t].append(s)
    basis = [tuple(int(i == j) for i in range(n)) for j in range(n)]
    beta0 = (r-1, r) + tuple(2*r-j//2 for j in range(2, n-1)) + (r,)
    delta = (2,4,6,5,4,3,2,1) + (0,)*(n-9) + (3,)
    gamma = (1,2,3,2,2,1,1,0) + (0,)*(n-9) + (1,)
    I = {0, n-1} | set(range(3, n-1, 2))

    def cartan(x):
        return tuple(2*x[j]-sum(x[t] for t in adj[j]) for j in range(n))

    def pair(x, y):
        return sum(v*w for v, w in zip(x, cartan(y)))

    def add(x, y, scale=1):
        return tuple(v+scale*w for v, w in zip(x, y))

    def refl(root, x):
        return add(x, root, -pair(root, x))

    gd = add(gamma, delta)

    def T(x):
        return refl(gamma, refl(gd, x))

    def N(x):
        return add(T(x), x, -1)

    # Uniform positive-real-root witness for gamma and gamma+delta, supported
    # in the same affine parabolic.  A greedy path is a certificate, not a
    # criterion claiming that norm two by itself implies a real root.
    witnesses = []
    for root in (gamma, gd):
        x, word = root, []
        while sum(x) > 1:
            p = cartan(x)
            s = next(j for j in range(n) if p[j] > 0)
            y = refl(basis[s], x)
            assert min(y) >= 0 and sum(y) < sum(x)
            assert s in set(range(8)) | {n-1}
            word.append("branch" if s == n-1 else s)
            x = y
        assert x in basis and pair(root, root) == 2
        witnesses.append({"lowering_word": word,
                          "ends_at": "branch" if x == basis[n-1] else basis.index(x)})
    assert pair(delta, delta) == pair(gamma, delta) == 0
    assert pair(beta0, gamma) == -r and pair(beta0, delta) == -a
    for e in basis:
        assert T(e) == add(add(e, gamma, pair(e,delta)), delta,
                           -pair(e,gamma)-pair(e,delta))
        assert N(N(N(e))) == (0,)*n
    v = tuple((a+r)*d-a*g for d, g in zip(delta,gamma))
    u = tuple(2*a*d for d in delta)
    assert N(beta0) == v and N(v) == u and N(u) == (0,)*n
    assert min(beta0) > 0 and min(v) >= 0 and min(u) >= 0
    # Independent D-coordinate realization, doubled to avoid fractions.
    after0 = refl(basis[0], beta0)
    finite = [-after0[0]]*(n-1)
    finite[0] += 2*after0[1]-2*after0[n-1]
    finite[1] += 2*after0[1]+2*after0[n-1]
    for j in range(2, n-1):
        finite[j] += 2*after0[j]
        finite[j-1] -= 2*after0[j]
    assert after0[0] == 1
    assert finite == [-1 if i%2 == 1 else 1 for i in range(1, n)]
    assert sum(c > 0 for c in finite) == 2*r  # allowed even sign changes
    matching = [(0,1),(2,n-1)] + [(j,j+1) for j in range(3,n-2,2)]
    assert len(matching) == 2*r and len(set(sum((list(e) for e in matching), []))) == 4*r
    assert len(I) == n-len(matching) == 2*r+1
    assert not any(s in I and t in I for s,t in edges)
    m0,m1,m2 = map(cartan,(beta0,v,u))
    assert all(m0[j]>0 and m1[j]>=0 and m2[j]>=0 for j in I)
    assert all(m0[j]<0 and m1[j]<=0 and m2[j]<=0 for j in range(n) if j not in I)
    for s,t in edges:
        if s in I or t in I:
            assert all(m[s]+m[t] <= 0 for m in (m0,m1,m2))
    rows=[]
    for k in (0,1,2,10,10**6):
        root=tuple(b-a*k*g+(a*k*k+r*k)*d for b,g,d in zip(beta0,gamma,delta))
        assert root == tuple(b+k*vj+k*(k-1)//2*uj for b,vj,uj in zip(beta0,v,u))
        assert min(root)>0 and pair(root,root)==2
        m=cartan(root)
        assert {j for j in range(n) if m[j]>0} == I
        columns=[refl(root,e) for e in basis]
        assert {j for j,col in enumerate(columns) if max(col)<=0} == I
        for s,t in edges:
            if s in I or t in I:
                assert min(add(columns[s],columns[t])) >= 0
        assert all(columns[0][j] != basis[0][j] for j in range(n))
        assert root[0] == r-1+4*k+(4*r-8)*k*k
        rows.append({"k":k,"height":sum(root),"coordinate_0":root[0]})
    return {"r":r,"rank":n,"alpha":len(I),"sample_checks":rows},witnesses


def main():
    data={"symbolic":symbolic(), "scope":"All-rank proof in referee.txt; matrix rows are finite counterchecks.",
          "independence":"No construction, targeted, or FC enumeration code is imported.","rows":[]}
    for r in range(3,41):
        row,witnesses=check_rank(r)
        if r==3:data["uniform_affine_root_witnesses"]=witnesses
        else:assert witnesses==data["uniform_affine_root_witnesses"]
        data["rows"].append(row)
    source=ROOT/"results/exceptional-leading.tex"
    data["audited_source_sha256"]=sha256(source.read_bytes()).hexdigest()
    (HERE/"referee-certificate.json").write_text(json.dumps(data,indent=2)+"\n")
    print("Symbolic identities and independent exact matrices passed: r=3..40, k=0,1,2,10,10^6.")
    print("Affinely supported real-root witnesses are identical in all checked ranks.")
    print("Source SHA256:",data["audited_source_sha256"])


if __name__=="__main__":main()
