#!/usr/bin/env python3
"""Exact check of mu(x6,w6) in the D6 parabolic of E7 (manuscript numbering), via the
KL recursion restricted to the Bruhat interval [e,w]. Also checks that the E7 word
132543621324356 equals the image of Gern's w6 = [2,0][4,0][6,0][5,4][6,6] under the
label map 1,2,3,4,5,6 -> 1,6,2,3,4,5, and that Gern's signed permutations match."""
import sys
from functools import lru_cache
sys.setrecursionlimit(10000)

# ---------- Gern's D6 as signed permutations (one-line notation, positions 1..6)
nD = 6
def gern_rmul(w, i):
    w = list(w)
    if i == 1:
        w[0], w[1] = -w[1], -w[0]
    else:
        w[i-2], w[i-1] = w[i-1], w[i-2]
    return tuple(w)
def gern_word_to_perm(word):
    w = tuple(range(1, nD+1))
    for i in word:
        w = gern_rmul(w, i)
    return w
def interval(j, i):  # Gern [j,i] with j>=i>=2 ... we only need [j,0] = (s1...sj)^{-1} and [5,4], [6,6]
    pass
# w6 = [2,0][4,0][6,0][5,4][6,6]; [j,0] = s_j ... s_2 s_1 ; [5,4] = ([4,5])^{-1} = s5 s4 ; [6,6] = s6
w6_word = [2,1] + [4,3,2,1] + [6,5,4,3,2,1] + [5,4] + [6]
x6_word = [1,2,4,6]
print("Gern w6 word:", w6_word, "length", len(w6_word))
print("Gern w6 signed perm:", gern_word_to_perm(w6_word))
print("Gern x6 signed perm:", gern_word_to_perm(x6_word))
# descent sets via Prop 2.2.4: R(w) = {s_i : w(i-1) > w(i)}, w(0) = -w(2)
def gern_R(w):
    ext = lambda k: -w[1] if k == 0 else w[k-1]
    return sorted(i for i in range(1, nD+1) if ext(i-1) > ext(i))
def gern_inv(w):
    inv = [0]*nD
    for pos, val in enumerate(w, start=1):
        inv[abs(val)-1] = pos if val > 0 else -pos
    return tuple(inv)
W6 = gern_word_to_perm(w6_word); X6 = gern_word_to_perm(x6_word)
print("R(w6)=", gern_R(W6), "L(w6)=", gern_R(gern_inv(W6)))
print("R(x6)=", gern_R(X6), "L(x6)=", gern_R(gern_inv(X6)))
# length via Prop 2.2.2
def gern_len(w):
    c = 0
    for i in range(nD):
        for j in range(i+1, nD):
            if w[i] > w[j]: c += 1
            if -w[i] > w[j]: c += 1
    return c
print("len(w6)=", gern_len(W6), "len(x6)=", gern_len(X6))

# ---------- E7 geometric representation (manuscript numbering), compare words
n = 7
adj = {i: set() for i in range(n)}
for i in range(n-2):
    adj[i].add(i+1); adj[i+1].add(i)
adj[n-1].add(2); adj[2].add(n-1)
def rmul(w, i):  # w: tuple of n columns (each tuple of n ints)
    cols = list(w)
    ai = cols[i]
    cols[i] = tuple(-c for c in ai)
    for j in adj[i]:
        cols[j] = tuple(a+b for a, b in zip(w[j], ai))
    return tuple(cols)
def ident():
    return tuple(tuple(1 if i == j else 0 for i in range(n)) for j in range(n))
def word_to_elt(word):
    w = ident()
    for s in word: w = rmul(w, s)
    return w
def neg(col): return any(c < 0 for c in col) and all(c <= 0 for c in col)
def R(w): return {j for j in range(n) if neg(w[j])}
def redword(w):
    word = []
    while True:
        r = R(w)
        if not r: break
        s = min(r); w = rmul(w, s); word.append(s)
    return word[::-1]
def length(w): return len(redword(w))
gmap = {1: 1, 2: 6, 3: 2, 4: 3, 5: 4, 6: 5}
w6_E7 = word_to_elt([gmap[s] for s in w6_word])
x6_E7 = word_to_elt([gmap[s] for s in x6_word])
tab_w = word_to_elt([int(c) for c in "132543621324356"])
tab_x = word_to_elt([int(c) for c in "1356"])
print("mapped Gern w6 == table E7 terminal:", w6_E7 == tab_w, " mapped x6 == table x:", x6_E7 == tab_x)
print("table word length:", length(tab_w))

# ---------- KL polynomials on the interval [e, w] in E7 (equivalently D6 parabolic)
# Polynomials as tuples of ints (coefficient list). Use the standard recursion with s in L(w):
# P_{x,w} = q^{1-c} P_{sx,sw} + q^c P_{x,sw} - sum_{z: sz<z, z<sw} mu(z,sw) q^{(l(w)-l(z))/2} P_{x,z},  c = 1 if sx<x else 0
def inverse(w):
    word = redword(w)
    return word_to_elt(word[::-1])
def L(w): return R(inverse(w))
def lmul(s, w):  # s w: apply s to each column
    cols = []
    for col in w:
        p = 2*col[s] - sum(col[j] for j in adj[s])
        c = list(col); c[s] -= p; cols.append(tuple(c))
    return tuple(cols)
# enumerate the lower Bruhat interval of tab_w via subwords is expensive; instead build downward closure:
# all elements below w = all subwords; equivalently BFS from w by "x -> x t for t with l(xt)<l(x)" does NOT give Bruhat interval.
# Use: Bruhat interval [e,w] = {x : x <= w}; test x<=w via subword property using a reduced word of w (greedy algorithm):
def bruhat_le(x, wword):
    # standard greedy: x <= w iff ... use the criterion: for w = s w', x <= w iff (sx<x ? sx <= w' : x <= w' ) -- Deodhar/lifting property variant:
    # if s in L(w): x <= w iff min(x, sx) <= sw.  (lifting property)
    if not wword: return x == ident()
    s = wword[0]; wprime = wword[1:]
    sx = lmul(s, x)
    if s in L(x):
        return bruhat_le(sx, wprime)
    else:
        return bruhat_le(x, wprime)
Wword = redword(tab_w)
Lmemo = {}
def Lc(w):
    if w not in Lmemo: Lmemo[w] = L(w)
    return Lmemo[w]
lenmemo = {}
def ln(w):
    if w not in lenmemo: lenmemo[w] = length(w)
    return lenmemo[w]
# Collect the interval: generate the subgroup elements by BFS from identity limited to length<=15 and test bruhat_le
elems = {ident()}
frontier = [ident()]
while frontier:
    nxt = []
    for w in frontier:
        if ln(w) >= len(Wword): continue
        for s in range(n):
            v = rmul(w, s)
            if v not in elems and ln(v) == ln(w)+1:
                elems.add(v); nxt.append(v)
    frontier = nxt
# restrict to D6 parabolic support automatically by Bruhat test
redw = {}
def rw(w):
    if w not in redw: redw[w] = redword(w)
    return redw[w]
below = [x for x in elems if bruhat_le(x, Wword)]
print("size of interval [e,w]:", len(below))
belowset = set(below)
def padd(a, b):
    m = max(len(a), len(b)); return tuple((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(m))
def psub(a, b):
    m = max(len(a), len(b)); return tuple((a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0) for i in range(m))
def pshift(a, k): return tuple([0]*k + list(a))
def pscale(a, c): return tuple(c*v for v in a)
def trim(a):
    a = list(a)
    while a and a[-1] == 0: a.pop()
    return tuple(a)
ZERO = (); ONE = (1,)
@lru_cache(maxsize=None)
def P(x, w):
    # requires x,w in the ambient group; returns KL polynomial P_{x,w}
    if not bruhat_le(x, rw(w)): return ZERO
    if x == w: return ONE
    s = min(Lc(w))
    sw = lmul(s, w)
    sx = lmul(s, x)
    c = 1 if s in Lc(x) else 0
    res = padd(pshift(P(sx, sw), 1-c), pshift(P(x, sw), c))
    lw = ln(w)
    for z in interval_below(sw):
        if s in Lc(z) and z != sw:
            m = mu(z, sw)
            if m:
                res = psub(res, pscale(pshift(P(x, z), (lw - ln(z))//2), m))
    return trim(res)
@lru_cache(maxsize=None)
def interval_below(w):
    ww = rw(w)
    return tuple(x for x in below if ln(x) <= len(ww) and bruhat_le(x, ww))
def mu(x, w):
    d = ln(w) - ln(x)
    if d <= 0 or d % 2 == 0: return 0
    p = P(x, w)
    k = (d-1)//2
    return p[k] if k < len(p) else 0
print("P_{x,w} =", P(tab_x, tab_w))
print("mu(x,w) =", mu(tab_x, tab_w), " l(w)-l(x) =", ln(tab_w)-ln(tab_x))
# sanity: P_{e,w}
print("P_{e,w} =", P(ident(), tab_w))
# also check descents of table terminal
print("L(w)=", sorted(Lc(tab_w)), "R(w)=", sorted(R(tab_w)), "L(x)=", sorted(Lc(tab_x)), "R(x)=", sorted(R(tab_x)))
