import sys, time
sys.path.insert(0, '.')
from klfinite import *
t0 = time.time()
# D_n lengths of x_n, w_n via formula for n=4..12 (no enumeration)
def xn(n):
    w = tuple(range(1, n+1))
    for s in [0, 1] + list(range(3, n, 2)):  # Gern s1,s2,s4,...,s_n -> labels 0,1,3,5,...,n-1
        w = act('D', w, s)
    return w
def wn(n):
    return tuple(((-1)**(n//2)) if i == 1 else (i if i % 2 else -(n+2-i)) for i in range(1, n+1))
for n in range(4, 25, 2):
    print(n, "l(x)=", length('D', xn(n)), "l(w)=", length('D', wn(n)), "formula", 3*n*n//8 + n//4, "neg entries", sum(1 for a in wn(n) if a<0))
G6 = Group('D', 6)
print("D6 order", len(G6.elems), time.time()-t0)
x, w = (-1,-2,4,3,6,5), (-1,-6,3,-4,5,-2)
assert x == xn(6) and w == wn(6)
I = G6.interval(x, w)
print("interval size", len(I), "lower ideal size", len(G6.ideal(w)))
rv = [0]*12
for z in I: rv[G6.len[z]-G6.len[x]] += 1
print("rank vector", rv)
P = G6.P(x, w); print("P(x6,w6) =", P, "mu =", G6.mu(x, w), time.time()-t0)
# all right descents
for s in G6.rdesc(w):
    G6._kl = {}
    # force recurrence with this descent: temporarily monkeypatch
    orig = G6.rdesc
    G6.rdesc = lambda v, s0=s, o=orig: ([s0] + [t for t in o(v) if t != s0]) if v == w else o(v)
    print(" descent", s, "->", G6.P(x, w))
    G6.rdesc = orig
print("reciprocity:", G6.check_reciprocity(x, w), time.time()-t0)
atoms = [z for z in I if G6.len[z] == 5]
print("atoms:", len(atoms))
a, b = (-4,-2,1,3,6,5), (-1,-3,4,2,6,5)
assert a in atoms and b in atoms
common = [z for z in I if G6.leq(a, z) and G6.leq(b, z)]
minimal = [z for z in common if not any(t != z and G6.leq(t, z) for t in common)]
print("common upper bounds", len(common), "minimal:", minimal, [G6.len[z]-4 for z in minimal])
print("c,d incomparable:", not G6.leq(minimal[0], minimal[1]) and not G6.leq(minimal[1], minimal[0]))
# P_{e,w6} equals P_{x6,w6}? (Lemma 2.3.9)
print("P(e,w6) =", G6.P(tuple(range(1,7)), w), time.time()-t0)
