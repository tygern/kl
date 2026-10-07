"""Independently verify the shape of Gern's I6 using direct subwords.

No SparseCoxeter import, lifting-property comparison, or KL computation.
Run from the project root: python3 research/verify_i6_subword.py
"""

from functools import lru_cache
import json


N = 6
IDENTITY = tuple(range(1, N + 1))
BOTTOM = (-1, -2, 4, 3, 6, 5)
TOP = (-1, -6, 3, -4, 5, -2)
ATOM_A = (-4, -2, 1, 3, 6, 5)
ATOM_B = (-1, -3, 4, 2, 6, 5)
EXPECTED_MINIMA = {
    (-4, -3, 1, 2, 6, 5),
    (-2, -4, 3, 1, 6, 5),
}


@lru_cache(None)
def length(v):
    """Type-D length by signed inversion counts."""
    return sum(v[i] > v[j] for i in range(N) for j in range(i + 1, N)) + sum(
        -v[i] > v[j] for i in range(N) for j in range(i + 1, N)
    )


def right(v, s):
    """Gern's simple generators, with zero-based labels."""
    u = list(v)
    if s == 0:
        u[0], u[1] = -u[1], -u[0]
    else:
        u[s - 1], u[s] = u[s], u[s - 1]
    return tuple(u)


@lru_cache(None)
def word(v):
    """Obtain one reduced word by removing a right descent."""
    if v == IDENTITY:
        return ()
    for s in range(N):
        z = right(v, s)
        if length(z) < length(v):
            return word(z) + (s,)
    raise AssertionError(v)


@lru_cache(None)
def lower(v):
    """Generate the lower ideal by all subwords, without order comparisons."""
    elements = {IDENTITY}
    for s in word(v):
        elements |= {right(z, s) for z in elements}
    return frozenset(elements)


def main():
    interval = {z for z in lower(TOP) if BOTTOM in lower(z)}
    atoms = [z for z in interval if length(z) == length(BOTTOM) + 1]
    common = [z for z in interval if ATOM_A in lower(z) and ATOM_B in lower(z)]
    minimal = [
        z for z in common if not any(t != z and t in lower(z) for t in common)
    ]
    assert len(interval) == 1676
    assert length(TOP) - length(BOTTOM) == 11
    assert len(atoms) == 12
    assert ATOM_A in atoms and ATOM_B in atoms
    assert len(common) == 720
    assert set(minimal) == EXPECTED_MINIMA
    print(json.dumps({
        "interval_size": len(interval),
        "rank": length(TOP) - length(BOTTOM),
        "atoms": len(atoms),
        "all_common_bounds": len(common),
        "minima": sorted(minimal),
        "relative_ranks": [length(z) - length(BOTTOM) for z in sorted(minimal)],
        "method": "Independent direct subword membership; no lifting-property comparisons or KL computation",
        "status": "All assertions passed",
    }, separators=(",", ":")))


if __name__ == "__main__":
    main()
