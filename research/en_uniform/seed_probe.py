"""Exact integral lowering of a symbolic rank-uniform root candidate.

This script discovers witnesses; finite successful ranks are not an all-rank
proof. No group or KL interval is enumerated.
"""
import json
from pathlib import Path


def data(r):
    n = 4*r+1
    adj = [[] for _ in range(n)]
    for a, b in [(i, i+1) for i in range(n-2)] + [(2, n-1)]:
        adj[a].append(b)
        adj[b].append(a)
    beta = [r-1, r, 2*r-1, 2*r-1]
    for t in range(2*r-2, 0, -1):
        beta.extend([t, t])
    beta.append(r)
    assert len(beta) == n
    return adj, beta


def lower(r):
    adj, beta = data(r)
    start = beta[:]
    steps = []
    while sum(beta) > 1:
        choices = [(i, 2*c-sum(beta[j] for j in adj[i]))
                   for i, c in enumerate(beta)]
        choices = [(i, m) for i, m in choices if 0 < m <= beta[i]]
        if not choices:
            return dict(r=r, rank=len(beta), success=False, stuck=beta,
                        steps=steps)
        i, m = choices[0]
        beta[i] -= m
        steps.append(i)
    return dict(r=r, rank=len(beta), success=True, seed=start,
                simple_root=beta.index(1), steps=steps)


if __name__ == '__main__':
    records = [lower(r) for r in range(3, 31)]
    path = Path(__file__).with_name('seed-probe.json')
    path.write_text(json.dumps(records, indent=2)+'\n')
    for row in records:
        print(row['r'], row['rank'], row['success'], len(row['steps']))
