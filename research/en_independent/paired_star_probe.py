"""Bounded simultaneous-star searches for the unresolved affine reflection pair.

A discovered descent mismatch/incomparability is a proof witness for mu=0;
a finite closed component without such a witness is not a mu computation.
Uses the saved research exact-matrix engine; independently audited elsewhere.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'research/broad_cells'))
from targeted import en

def main():
    g = en(9)
    seed = json.loads((ROOT / 'research/en_families/affine_reflection_family.json').read_text())
    wbase = g.elt(seed['base_reflection_word'])
    edges = g.edges
    def star(w, side, s, t):
        a = g.inv(w) if side else w
        if len(set(g.desc(a)) & {s,t}) != 1:
            return None
        moves = [g.right(a, u) for u in (s,t)
                 if len(set(g.desc(g.right(a,u))) & {s,t}) == 1]
        assert len(moves) == 1
        return g.inv(moves[0]) if side else moves[0]
    summaries = []
    for extra in (None,0,1):
        w = wbase if extra is None else g.right(wbase,extra)
        # Bare pair includes odd-gap bottom of length4; extensions use max I.
        x = g.elt(g.desc(w))
        start = (x,w)
        todo = [start]
        index = {start:0}
        parent = [(None,None)]
        witness = None
        for i,(a,b) in enumerate(todo):
            la,lb = g.length(a),g.length(b)
            low,high = (a,b) if la < lb else (b,a)
            gap = abs(lb-la)
            if gap == 1 and g.leq(low,high):
                witness = {'kind':'cover_mu_one','state':i}
                break
            if not g.leq(low,high):
                witness = {'kind':'incomparable_mu_zero','state':i}
                break
            if gap > 1 and (not set(g.desc(high)) <= set(g.desc(low)) or
                            not set(g.desc(g.inv(high))) <= set(g.desc(g.inv(low)))):
                witness = {'kind':'descent_mismatch_mu_zero','state':i}
                break
            for side in (0,1):
                for s,t in edges:
                    aa,bb = star(a,side,s,t),star(b,side,s,t)
                    if aa is None or bb is None:
                        continue
                    pair = (aa,bb)
                    if pair not in index:
                        index[pair] = len(todo)
                        todo.append(pair)
                        parent.append((i,(side,s,t)))
            if len(todo) >= 100000:
                break
        row = {'extra_generator':extra,'start_lower_word':g.word(x),'start_upper_word':g.word(w),
               'states_discovered':len(todo),'states_processed':i+1,
               'closed':witness is None and i+1==len(todo),'witness':witness}
        if witness:
            steps=[];j=witness['state']
            while parent[j][0] is not None:
                p,move=parent[j];steps.append(move);j=p
            row['witness_steps']=list(reversed(steps))
            a,b=todo[witness['state']]
            row['witness_pair']=[g.word(a),g.word(b)]
        summaries.append(row)
        print(json.dumps(row),flush=True)
        g.right.cache_clear();g.desc.cache_clear();g.word.cache_clear();g.inv.cache_clear()
    (ROOT / 'research/en_independent/paired_star_probe.json').write_text(json.dumps(summaries,indent=2)+'\n')

if __name__=='__main__':
    main()
