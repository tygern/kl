"""Compare the independent finite terminal certificates and known invariants.

python3 research/verify_exceptional.py             # audit saved certificates
python3 research/verify_exceptional.py --recompute # rerun all exact engines

Four snapshots are compared: the root-index enumerations of E6 and E7
(research/broad_exceptional/e6_bad.json and e7_bad.json, both written by
enumerate_bad.cpp) and the integer-matrix enumerations (results/
e6-independent-certificate.json by verify_e6_independent.py and results/
e7-independent-certificate.json by verify_e7_matrices.cpp).  The two
programs of each pair use different representations and share no code.
"""
import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def run(args):
    return subprocess.run(args, cwd=ROOT, check=True, text=True, capture_output=True).stdout

def read(path):
    return json.loads((ROOT / path).read_text())

def polynomial(exponents):
    result = [1]
    for exponent in exponents:
        new = [0] * (len(result) + exponent)
        for i, c in enumerate(result):
            for j in range(exponent + 1):
                new[i + j] += c
        result = new
    return result

def canonical_rootindex(rows):
    """Rows of enumerate_bad.cpp: word, length, Rmask, Lmask, bottoms[word, length, rank]."""
    return sorted((tuple(row['word']), row['length'],
                   tuple((tuple(x['word']), x['length'], x['rank']) for x in row['bottoms']))
                  for row in rows)

def canonical_e6_matrix(rows):
    """Rows of verify_e6_independent.py: word, length, candidates[word, length, interval_rank]."""
    return sorted((tuple(row['word']), row['length'],
                   tuple((tuple(x['word']), x['length'], x['interval_rank']) for x in row['candidates']))
                  for row in rows)

def canonical_e7_matrix(rows):
    """Rows of verify_e7_matrices.cpp: word, length, bottoms[word, length, rank]."""
    return canonical_rootindex(rows)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--recompute', action='store_true')
    args = parser.parse_args()
    if args.recompute:
        compiler = next((shutil.which(c) for c in ('c++', 'clang++', 'g++') if shutil.which(c)), None)
        if compiler is None:
            raise SystemExit('A C++17 compiler is required for the E6/E7 reproductions.')
        run([sys.executable, str(ROOT / 'research/verify_e6_independent.py')])
        with tempfile.TemporaryDirectory(prefix='kl-exceptional-') as folder:
            origin = str(Path(folder) / 'root-ids')
            matrix = str(Path(folder) / 'matrices')
            run([compiler, '-std=c++17', '-O3', str(ROOT / 'research/broad_exceptional/enumerate_bad.cpp'), '-o', origin])
            run([compiler, '-std=c++17', '-O3', str(ROOT / 'research/verify_e7_matrices.cpp'), '-o', matrix])
            (ROOT / 'research/broad_exceptional/e6_bad.json').write_text(run([origin, '6']))
            (ROOT / 'research/broad_exceptional/e7_bad.json').write_text(run([origin, '7']))
            (ROOT / 'results/e7-independent-certificate.json').write_text(run([matrix]))
    a6 = read('research/broad_exceptional/e6_bad.json')
    b6 = read('results/e6-independent-certificate.json')
    a7 = read('research/broad_exceptional/e7_bad.json')
    b7 = read('results/e7-independent-certificate.json')
    assert canonical_rootindex(a6['bad']) == canonical_e6_matrix(b6['non_fc_weak_bad_terminals'])
    assert canonical_rootindex(a7['bad']) == canonical_e7_matrix(b7['bad'])
    # Poincare polynomials of E6 and E7 from their exponents, element counts,
    # fully commutative counts (Stembridge 1998) and commuting terminal counts.
    e6 = polynomial([1, 4, 5, 7, 8, 11])
    e7 = polynomial([1, 5, 7, 9, 11, 13, 17])
    assert sum(e6) == 51840 and sum(e7) == 2903040
    assert a6['length_distribution'] == [b6['length_distribution'][str(i)] for i in range(len(e6))] == e6
    assert a7['length_distribution'] == b7['length_distribution'] == e7
    assert a6['order'] == b6['element_count'] == 51840
    assert a7['order'] == b7['order'] == 2903040
    assert a6['fc_count'] == b6['fc_count'] == 662
    assert a7['fc_count'] == b7['fc_count'] == 2670
    assert a6['commuting_weak'] == b6['commuting_terminals'] == 22
    assert a7['commuting_weak'] == b7['commuting_terminals'] == 36
    # Star-operation closure checks: the matrix programs count both sides.
    assert b6['fc_star_checks'] == 2 * a6['fc_right_star_checks'] == 3620
    assert b7['fc_star_checks'] == 2 * a7['fc_right_star_checks'] == 16776
    # Independent signed permutation check of the odd terminal, without a KL
    # evaluator or imports from either exceptional enumeration.
    mapping = {1: 0, 6: 1, 2: 2, 3: 3, 4: 4, 5: 5}
    def signed(word):
        values = list(range(1, 7))
        for label in word:
            s = mapping[label]
            if s == 0:
                values[0], values[1] = -values[1], -values[0]
            else:
                values[s-1], values[s] = values[s], values[s-1]
        return values
    odd = [row for row in b7['bad'] if row['length'] == 15][0]
    assert signed(odd['word']) == [-1, -6, 3, -4, 5, -2]
    assert signed(odd['bottoms'][0]['word']) == [-1, -2, 4, 3, 6, 5]
    paths = [
        'research/verify_e6_independent.py', 'research/verify_e7_matrices.cpp',
        'research/broad_exceptional/exact_e.py', 'research/broad_exceptional/enumerate_bad.cpp',
        'results/e6-independent-certificate.json', 'results/e7-independent-certificate.json',
        'research/broad_exceptional/e6_bad.json', 'research/broad_exceptional/e7_bad.json',
    ]
    report = {'status': 'passed', 'independent_terminal_tables_agree': True,
              'poincare_histograms_match': True, 'd6_signed_pair_matches': True,
              'snapshot_generators': {'research/broad_exceptional/e6_bad.json': 'enumerate_bad.cpp 6',
                                      'research/broad_exceptional/e7_bad.json': 'enumerate_bad.cpp 7',
                                      'results/e6-independent-certificate.json': 'verify_e6_independent.py',
                                      'results/e7-independent-certificate.json': 'verify_e7_matrices.cpp'},
              'files': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths},
              'priority': 'Not established by this computational audit.'}
    (ROOT / 'results/exceptional-audit.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))

if __name__ == '__main__':
    main()
