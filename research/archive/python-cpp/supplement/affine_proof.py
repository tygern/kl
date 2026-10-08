"""Proof-only entry point using the unchanged affine referee arithmetic.

Unlike verify_families.main(), this does not audit exploratory secondary families.
Run from the supplement runner; its working tree has the original source layout.
"""
import json
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'research/en_affine_referee'))
import verify_families as q


def main():
    assert all(q.pair(q.DELTA, e) == 0 for e in q.E)
    assert q.GAMMA in q.FINITE_ROOTS
    assert q.mm(q.reflection(q.GAMMA), q.reflection(q.add(q.GAMMA, q.DELTA))) == q.translate(q.GAMMA)
    data = json.loads((ROOT/'research/en_families/e9_full_support_eligible.json').read_text())
    word = next(row['word'] for row in data if row['length'] == 27)
    a = q.word_matrix(word)
    assert q.mm(a,a) == q.E
    d = tuple(q.pair(e,q.add(q.GAMMA,q.mv(a,q.GAMMA),-1)) for e in q.E)
    assert d == (2,-1,1,-1,1,-1,1,-1,-1)
    assert q.mm(q.mm(q.translate(q.GAMMA),a),q.translate(q.GAMMA,-1)) == q.shifted(a,d,1)
    result = q.audit_family(word,d,(27,92),[1,3,5,7,8],'proved odd-length affine conjugates')
    assert result['complete_FC_mask_candidates'] == 1
    assert result['eligible_FC_bottoms_for_all_k'][0]['word'] == [1,3,5,7,8]
    (q.HERE/'proved-affine-certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'affine_length':[27,92], 'finite_roots':len(q.FINITE_ROOTS),
                      'terminal_for_all_k':result['terminal_for_all_k']},indent=2))

if __name__ == '__main__': main()
