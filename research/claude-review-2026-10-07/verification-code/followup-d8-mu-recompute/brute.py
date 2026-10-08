"""Independent brute-force KL polynomials for the whole group W(D_n), n small.
Uses Gern's Prop 1.2.6 literally (LEFT descents, left multiplication), subword
Bruhat criterion, and plain dict polynomials. Exact Python integers.
Output: canonical words (min left-descent peeling) and P_{x,y} for all x<=y.
"""
import sys
from itertools import permutations, product
from functools import lru_cache

n = int(sys.argv[1])
out = sys.argv[2]

def rmul(w, j):  # Gern s_j on the right, j=1..n
    w = list(w)
    if j == 1: w[0], w[1] = -w[1], -w[0]
    else: w[j-2], w[j-1] = w[j-1], w[j-2]
    return tuple(w)

def lmul(w, j):  # Gern s_j on the left: act on values
    def f(v):
        a, s = abs(v), (1 if v > 0 else -1)
        if j == 1:
            if a == 1: return -2*s
            if a == 2: return -1*s
        else:
            if a == j-1: return j*s
            if a == j: return (j-1)*s
        return v
    return tuple(f(v) for v in w)

def length(w):
    return sum(1 for i in range(n) for j in range(i+1, n) if w[i] > w[j]) + \
           sum(1 for i in range(n) for j in range(i+1, n) if w[i]+w[j] < 0)

# enumerate group
G = []
for p in permutations(range(1, n+1)):
    for signs in product((1, -1), repeat=n):
        if sum(1 for s in signs if s < 0) % 2 == 0:
            G.append(tuple(a*s for a, s in zip(p, signs)))
G.sort(key=length)
L = {w: length(w) for w in G}
e = tuple(range(1, n+1))
assert len(G) == 2**(n-1) * __import__('math').factorial(n)

def ldes(w): return [j for j in range(1, n+1) if L[lmul(w, j)] < L[w]]
def rdes(w): return [j for j in range(1, n+1) if L[rmul(w, j)] < L[w]]
LD = {w: ldes(w) for w in G}

# canonical word by min left descent
canon = {e: ''}
for w in G:
    if w == e: continue
    j = min(LD[w]); canon[w] = str(j) + canon[lmul(w, j)]

# reduced word (right-peeling) and subword Bruhat test
def redword(w):
    word = []
    while w != e:
        j = min(rdes(w)); word.append(j); w = rmul(w, j)
    return word[::-1]
RW = {w: redword(w) for w in G}

def is_subword_of(x, word):
    # does some subword of `word` multiply (right-to-left applied as right mult) to x?  Greedy from the left:
    # standard: x <= w iff x is a subword of a reduced word of w; test via DP over positions.
    states = {e}
    for j in word:
        new = set(states)
        for s in states:
            new.add(rmul(s, j))
        states = new
    return x in states
leq = {}
for w in G:
    word = RW[w]
    states = {e}
    for j in word:
        states |= {rmul(s, j) for s in states}
    for x in states:
        leq[(x, w)] = True

def padd(a, b):
    r = dict(a)
    for k, v in b.items(): r[k] = r.get(k, 0) + v
    return {k: v for k, v in r.items() if v}
def pshift(a, k, c=1): return {d+k: c*v for d, v in a.items()}

P = {}
mu = {}
def getP(x, w): return P.get((x, w), {})
for w in G:
    if w == e:
        P[(e, e)] = {0: 1}; continue
    s = min(LD[w]); sw = lmul(w, s)
    # mu-list of sw: z with s z < z, z < sw, mu(z,sw) != 0
    zs = [(z, mu[(z, sw)]) for z in G if (z, sw) in leq and z != sw and L[lmul(z, s)] < L[z] and mu.get((z, sw), 0)]
    for x in G:
        if (x, w) not in leq: continue
        c = 1 if s in LD[x] else 0
        sx = lmul(x, s)
        p = padd(pshift(getP(sx, sw), 1-c), pshift(getP(x, sw), c))
        for z, m in zs:
            if (x, z) in leq:
                p = padd(p, pshift(getP(x, z), (L[w]-L[z])//2, -m))
        P[(x, w)] = p
        d = L[w]-L[x]
        if d % 2 == 1:
            mu[(x, w)] = p.get((d-1)//2, 0)
        else:
            mu[(x, w)] = 0

def pstr(p):
    if not p: return '0'
    s = ''
    for i in sorted(p):
        a = p[i]
        s += ('' if not s else (' - ' if a < 0 else ' + ')) + ('' if (abs(a) == 1 and i > 0) else str(abs(a))) + ('' if i == 0 else ('q' if i == 1 else 'q^%d' % i))
    return s
with open(out, 'w') as f:
    for (x, w), p in P.items():
        f.write('%s\t%s\t%s\n' % (canon[x], canon[w], pstr(p)))
print('group', len(G), 'pairs', len(P))
