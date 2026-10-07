"""Certify all surviving rank-vector candidates extend to lower-ideal isos.

Every mask of a fixed word is generated simultaneously in both systems.
Bijectivity of the resulting relation and preservation of all element lengths
give an order isomorphism by the subword property, in both directions.
"""
from coxeter_surgery import Coxeter, d_edges
import json
from pathlib import Path


def main():
    rows = [json.loads(line) for line in Path(__file__).with_name(
        'lift-all-bottoms-results.jsonl').read_text().splitlines()]
    for row in rows:
        if 'bottom' not in row: continue
        source = Coxeter(6, d_edges(6))
        edges = d_edges(6)
        i,j,_ = edges[row['edge']]
        edges[row['edge']] = (i,j,2)
        target = Coxeter(6,edges)
        pairs = {(source.e, target.e)}
        for s in row['top']:
            pairs |= {(source.right(a,s), target.right(b,s)) for a,b in pairs}
        source_values = {a for a,b in pairs}
        target_values = {b for a,b in pairs}
        relation_bijective = len(source_values) == len(target_values) == len(pairs)
        length_preserved = all(source.length(a) == target.length(b) for a,b in pairs)
        print(json.dumps({'edge':row['edge'],'top_word':row['top'],
                          'pair_count':len(pairs),'source_count':len(source_values),
                          'target_count':len(target_values),
                          'relation_bijective':relation_bijective,
                          'all_lengths_preserved':length_preserved}),flush=True)
        assert relation_bijective and length_preserved


if __name__ == '__main__': main()
