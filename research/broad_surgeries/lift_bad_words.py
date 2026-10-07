"""Lift all reduced expressions of a D bad element across one inflated bond."""
from coxeter_surgery import Coxeter, d_edges, bad_word
from functools import lru_cache
import json


def run(n, edge):
    source = Coxeter(n, d_edges(n))
    edges = d_edges(n)
    i,j,_ = edges[edge]
    edges[edge] = (i,j,2)
    target = Coxeter(n, edges)
    top = source.element(bad_word(n))
    botword = (0,1) + tuple(range(3,n,2))
    src_ranks = source.ranks(source.interval(source.element(botword),top))
    src_lower_ranks = source.ranks(source.lower(top))
    @lru_cache(None)
    def lifts(v):
        if v == source.e:
            return {target.e:()}
        result = {}
        for s in range(n):
            if not source.positive(v,s):
                for z,word in lifts(source.right(v,s)).items():
                    result[target.right(z,s)] = word + (s,)
        return result
    lifted = lifts(top)
    print(json.dumps({"n":n,"edge":edge,"lift_count":len(lifted)}),flush=True)
    for v,word in lifted.items():
        ranks = target.ranks(target.interval(target.element(botword),v))
        if ranks == src_ranks:
            lowerranks = target.ranks(target.lower(v))
            print(json.dumps({"n":n,"edge":edge,"candidate_word":word,
                              "ranks":ranks,"principal_ranks_equal":lowerranks==src_lower_ranks}),flush=True)


if __name__ == '__main__':
    for edge in range(5):
        run(6,edge)
