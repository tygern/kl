"""Independent affine-D4 witness check via subwords and R reciprocity.

Uses faithful affine signed permutations, no project imports, no KL descent
recurrence or lifting comparison. All lower ideals are direct subword sets.
python3 research/verify_affine_d4_r.py
"""
from functools import cache
from itertools import product
from collections import Counter
from pathlib import Path
import json

E = ((1, 2, 3, 4), (0, 0, 0, 0))
QX = (0, 1, 3, 4)
QW = QX + (2,) + QX

def simple(w, s):
    p, t = map(list, w)
    if s == 0:
        p[0], p[1] = -p[1], -p[0]
    elif s < 4:
        p[s-1], p[s] = p[s], p[s-1]
    else:
        for a in p[-2:]:
            t[abs(a)-1] += 1 if a > 0 else -1
        p[-2], p[-1] = -p[-1], -p[-2]
    return tuple(p), tuple(t)

def evaluate(word):
    w = E
    for s in word:
        w = simple(w, s)
    return w

def subwords(word):
    for mask in product((0, 1), repeat=len(word)):
        yield tuple(s for s, b in zip(word, mask) if b)

def add(a, b, shift=0, scale=1):
    result = list(a) + [0]*max(0, len(b)+shift-len(a))
    for k, c in enumerate(b):
        result[k+shift] += scale*c
    while result and result[-1] == 0:
        result.pop()
    return tuple(result)

def multiply(a, b):
    result = ()
    for k, c in enumerate(b):
        result = add(result, a, k, c)
    return result

def main():
    words = {}
    for word in subwords(QW):
        w = evaluate(word)
        if w not in words or len(word) < len(words[w]):
            words[w] = word
    lengths = {w: len(word) for w, word in words.items()}
    ideals = {w: frozenset(evaluate(word) for word in subwords(q)) for w, q in words.items()}
    def below(x, w):
        return x in ideals.get(w, ())
    def desc(w):
        return tuple(s for s in range(5) if lengths.get(simple(w,s), 100) < lengths[w])
    @cache
    def r(x, w):
        if not below(x, w):
            return ()
        if x == w:
            return (1,)
        s = desc(w)[0]
        xs, ws = simple(x, s), simple(w, s)
        if lengths.get(xs,100) < lengths[x]:
            return r(xs, ws)
        return add(add((), r(x,ws), scale=-1), add(r(x,ws), r(xs,ws)), shift=1)
    w = evaluate(QW)
    x = evaluate(QX)
    assert lengths[x] == 4 and lengths[w] == 9
    p = {w: (1,)}
    # q^d P(q^-1) - P(q) = sum_{x<z<=w} R_x,z(q) P_z,w(q).
    for y in sorted(words, key=lambda z: lengths[z], reverse=True):
        if y == w:
            continue
        d = lengths[w]-lengths[y]
        rhs = ()
        for z in words:
            if lengths[z] > lengths[y] and below(y,z):
                rhs = add(rhs, multiply(r(y,z), p[z]))
        bound = (d-1)//2
        value = tuple(-rhs[k] if k<len(rhs) else 0 for k in range(bound+1))
        while value and value[-1] == 0:
            value = value[:-1]
        p[y] = value
        dual = [0]*(d+1)
        for k, c in enumerate(value):
            dual[d-k] += c
        assert add(dual, value, scale=-1) == rhs, (words[y], value, rhs)
        assert value[0] == 1 and all(c>=0 for c in value)
    assert p[x] == (1,3,2), p[x]
    rank = dict(sorted(Counter(lengths[z]-lengths[x] for z in ideals[w] if below(x,z)).items()))
    report = {"method": "Affine signed actions, direct subwords, R-polynomial reciprocity",
              "simple_labels": "Forks0,1--center2--forks3,4; affine generator4",
              "bottom_word": QX, "top_word": QW, "bottom": x, "top": w,
              "lengths": [lengths[x], lengths[w]], "lower_ideal_size": len(words),
              "interval_rank_counts": rank, "polynomial": p[x], "mu": 2,
              "fc_bottom": "Four pairwise commuting generators, each once",
              "reciprocity_identities_checked": len(p)-1,
              "status": "passed"}
    output = Path(__file__).resolve().parents[1]/'results/affine-d4-independent.json'
    output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__ == '__main__':
    main()
