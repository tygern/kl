"""Check kl-transfer Section 5 E7/E8 coset example: a = w0(E_r) w0(J), J = {2..7} (Bourbaki), lengths and non-FC of a*x6, a*w6."""
import sys, time
sys.path.insert(0, '.')
from klfinite import Group
def cartan(r):
    # Bourbaki E_r: chain 1-3-4-5-6-7(-8), 2 attached to 4
    N = r; C = [[0]*N for _ in range(N)]
    for i in range(N): C[i][i] = 2
    edges = [(1,3),(3,4),(4,5),(5,6),(6,7),(2,4)] + ([(7,8)] if r == 8 else [])
    for a, b in edges: C[a-1][b-1] = C[b-1][a-1] = -1
    return C
def run(r, cap=3_000_000):
    t0 = time.time()
    C = cartan(r); N = r
    E = tuple(tuple(1 if i == j else 0 for i in range(N)) for j in range(N))
    def refl(s, v):
        coef = sum(C[s][j]*v[j] for j in range(N))
        return tuple(v[j] - (coef if j == s else 0) for j in range(N))
    def apply(w, v):
        out = [0]*N
        for j in range(N):
            if v[j]:
                for i in range(N): out[i] += v[j]*w[j][i]
        return tuple(out)
    def mul_s(w, s): return tuple(apply(w, refl(s, E[j])) for j in range(N))
    def neg(v): return all(c <= 0 for c in v) and any(c < 0 for c in v)
    def rdesc(w): return [s for s in range(N) if neg(w[s])]
    def length(w):
        n = 0
        while True:
            d = rdesc(w)
            if not d: return n
            w = mul_s(w, d[0]); n += 1
    w0 = tuple(tuple(-1 if i == j else 0 for i in range(N)) for j in range(N))
    L0 = length(w0); print(f"E{r}: l(w0) = {L0}")
    J = [1,2,3,4,5,6]  # labels 2..7 zero-based
    a = w0; k = 0
    while True:
        d = [s for s in rdesc(a) if s in J]
        if not d: break
        a = mul_s(a, d[0]); k += 1
    print(f"  stripped {k} J-descents; l(a) = {length(a)}; a in W^J: {not any(s in J for s in rdesc(a))}")
    # Gern D6 words for x6,w6 (code labels 0..5) -> E labels: code L -> Gern L+1 -> Bourbaki L+2 -> zero-based L+1
    G6 = Group('D', 6)
    x6, w6 = (-1,-2,4,3,6,5), (-1,-6,3,-4,5,-2)
    wx = [l+1 for l in G6.reduced_word(x6)]; ww = [l+1 for l in G6.reduced_word(w6)]
    def ev(w, word):
        for s in word: w = mul_s(w, s)
        return w
    ax, aw = ev(a, wx), ev(a, ww)
    print(f"  l(ax) = {length(ax)}, l(aw) = {length(aw)}")
    def noncomm(s, t): return C[s][t] != 0 and s != t
    def nonfc_witness(w):
        """Search prefixes u (right weak order) of w for two noncommuting right descents."""
        seen = {w}; stack = [w]
        while stack:
            u = stack.pop()
            d = rdesc(u)
            for i in range(len(d)):
                for j in range(i+1, len(d)):
                    if noncomm(d[i], d[j]): return True, len(seen)
            for s in d:
                v = mul_s(u, s)
                if v not in seen:
                    seen.add(v); stack.append(v)
            if len(seen) > cap: return None, len(seen)
        return False, len(seen)
    for name, el in [("a", a), ("ax", ax), ("aw", aw)]:
        res, n = nonfc_witness(el)
        print(f"  {name}: non-FC witness found = {res} (explored {n} prefixes), {time.time()-t0:.1f}s")
run(7)
run(8, cap=1_500_000)
