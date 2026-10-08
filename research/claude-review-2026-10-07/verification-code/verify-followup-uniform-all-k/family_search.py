import itertools, sys
from fc_check import adjacency_E, is_fc_reduced, length_by_matrix

# validate checker on enumerated max-length E13 word
w13 = [int(x) for x in "11 10 9 8 7 6 5 4 3 2 1 0 12 2 1 3 2 4 3 5 4 6 5 7 6 8 7 9 8 10 9 11 10 12 2 1 0 3 2 1 4 3 2 5 4 3 6 5 4 7 6 5 8 7 6 9 8 7 12 2 1 0 3 2 1 4 3 2 5 4 3 6 5 4 12 2 1 0 3 2 1 12 2 3 4 5 6 7 8 9 10 11".split()]
adj = adjacency_E(13)
print("E13 max word: FC&reduced =", is_fc_reduced(w13, adj), "matrix length =", length_by_matrix(w13, 13), "len =", len(w13))
# a non-reduced / non-FC word must fail
print("sanity: [0,1,0] FC?", is_fc_reduced([0,1,0], adj), " [0,0]?", is_fc_reduced([0,0], adj), " [0,1,2,0,1]?", is_fc_reduced([0,1,2,0,1], adj))

def run(a, b):
    return list(range(a, b+1))

# Search over interpretations of "blocks of width-2 runs for tops T with bottoms 12/012":
# a block is a concatenation of ascending runs with the given bottoms (in some order) and
# tops in {T, T-1, T-2} and an optional trailing p; we look for FC words of length 91 at n=13.
def build(n, blocktype):
    p = n-1
    A = list(range(n-2, 1, -1)) + [p]
    B1 = [1,2,3,0,1,2,p]
    tops = list(range(5, n-1, 2))
    word = A + B1
    nb = len(tops)
    for i, T in enumerate(tops):
        # last block is "12"; alternate backwards
        typ = "12" if (nb-1-i) % 2 == 0 else "012"
        word += blocktype(T, typ, p)
    word += run(0, n-2)
    return word

def make_blocktype(order12, tops12, order012, tops012, trailing):
    def bt(T, typ, p):
        if typ == "12":
            bottoms, tps = order12, tops12
        else:
            bottoms, tps = order012, tops012
        w = []
        for b, dt in zip(bottoms, tps):
            w += run(b, T-dt)
        if trailing: w.append(p)
        return w
    return bt

found = []
for order12 in itertools.permutations([1,2]):
    for tops12 in itertools.product([0,1,2], repeat=2):
        for order012 in itertools.permutations([0,1,2]):
            for tops012 in itertools.product([0,1,2], repeat=3):
                for trailing in (True, False):
                    bt = make_blocktype(order12, tops12, order012, tops012, trailing)
                    ok = True; lens = []
                    for n in (13, 17, 21):
                        w = build(n, bt)
                        if not is_fc_reduced(w, adjacency_E(n)):
                            ok = False; break
                        lens.append(len(w))
                    if ok:
                        found.append((order12, tops12, order012, tops012, trailing, lens))
print("interpretations that are FC for n=13,17,21:", len(found))
for f in found:
    print(f)
