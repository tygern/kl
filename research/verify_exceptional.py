"""Compare the independent finite terminal certificates and known invariants.

python3 research/verify_exceptional.py             # audit saved certificates
python3 research/verify_exceptional.py --recompute # rerun all exact engines
UV_CACHE_DIR=/private/tmp/kl-uv-cache uv run --offline --no-project \
  --python /opt/homebrew/bin/python3 research/verify_exceptional.py
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

def canonical(rows):
    def fields(row):
        bottoms = row.get('bottoms', row.get('xs', row.get('candidates')))
        return (tuple(row['word']), row['length'],
                tuple((tuple(x['word']), x['length'], x.get('rank', x.get('interval_rank')))
                      for x in bottoms))
    return sorted(fields(row) for row in rows)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--recompute', action='store_true')
    args = parser.parse_args()
    if args.recompute:
        # The two C++ engines use different representations and source files.
        compiler = shutil.which('c++')
        if compiler is None:
            raise SystemExit('A C++17 compiler is required for E7 reproduction.')
        run([sys.executable, str(ROOT / 'research/verify_e6_independent.py')])
        with tempfile.TemporaryDirectory(prefix='kl-exceptional-') as folder:
            origin = str(Path(folder) / 'root-ids')
            matrix = str(Path(folder) / 'matrices')
            run([compiler, '-std=c++17', '-O3', str(ROOT / 'research/broad_exceptional/enumerate_bad.cpp'), '-o', origin])
            run([compiler, '-std=c++17', '-O3', str(ROOT / 'research/verify_e7_matrices.cpp'), '-o', matrix])
            (ROOT / 'research/broad_exceptional/e7_bad.json').write_text(run([origin, '7']))
            (ROOT / 'results/e7-independent-certificate.json').write_text(run([matrix]))
    a6 = read('research/broad_exceptional/e6_bad.json')
    b6 = read('results/e6-independent-certificate.json')
    a7 = read('research/broad_exceptional/e7_bad.json')
    b7 = read('results/e7-independent-certificate.json')
    assert canonical(a6['bad']) == canonical(b6['non_fc_weak_bad_terminals'])
    assert canonical(a7['bad']) == canonical(b7['bad'])
    for rank, a, b, exponents, order, fc in (
        (6, a6, b6, [1, 4, 5, 7, 8, 11], 51840, 662),
        (7, a7, b7, [1, 5, 7, 9, 11, 13, 17], 2903040, 2670),
    ):
        expected = polynomial(exponents)
        def histogram(row):
            values = row['length_distribution']
            return values if isinstance(values, list) else [values[str(i)] for i in range(len(expected))]
        assert histogram(a) == histogram(b) == expected
        assert sum(expected) == order
        assert a['order'] == b.get('order', b.get('element_count')) == order
        assert a.get('fc_count', a.get('fc')) == b['fc_count'] == fc
        assert a['commuting_weak'] == b['commuting_terminals'] == (22 if rank == 6 else 36)
    assert b7['fc_star_checks'] == 2*a7['fc_right_star_checks'] == 16776
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
              'files': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths},
              'priority': 'Not established by this computational audit.'}
    (ROOT / 'results/exceptional-audit.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))

if __name__ == '__main__':
    main()
