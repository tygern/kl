"""Reproduce all finite checks and explicit isomorphism certificates.

Run from repository root: python3 computations/run_search.py
No third-party packages required; also works with `uv run`.
"""
import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path
from time import perf_counter

from coxeter import Coxeter, certify, isomorphism, poset_data, non_type_a_coset_certificate
from sparse_kl import SparseCoxeter, bad_d6
from verify_certificate import verify


def main():
    started = perf_counter()
    groups = {name: Coxeter(name[0], int(name[1:])) for name in ['A3', 'A4', 'B3', 'D4', 'B4']}
    result = dict(conventions=dict(
        permutations='one-line images; signed permutations act on {-n,...,-1,1,...,n}',
        multiplication='words multiply on the right; simple-generator labels start at 0',
        A='s_i swaps positions i and i+1, zero-based',
        B='s_0 negates first entry; s_i swaps entries i and i+1 for i>=1, one-based',
        D='s_0 sends first two entries (a,b) to (-b,-a); s_i swaps entries i and i+1 for i>=1, one-based',
        polynomial='coefficient arrays in ascending powers of q'),
        finite_checks={}, certificates=[], search_counts={})
    for name, group in groups.items():
        result['finite_checks'][name] = group.validate()
        sparse = SparseCoxeter(group.kind, group.rank)
        for w in range(group.order):
            wp = group.elements[w]
            assert sparse.length(wp) == group.length[w]
            for x in range(group.order):
                xp = group.elements[x]
                assert sparse.leq(xp, wp) == (x in group.lower[w])
                if x in group.lower[w]:
                    assert sparse.kl(xp, wp) == group.kl(x, w)
        result['finite_checks'][name]['checks'] += [
            'independent inversion-formula lengths for every element',
            'independent lifting-property Bruhat comparisons for every ordered pair',
            'second sparse KL implementation for every comparable pair']
        print(name, result['finite_checks'][name], flush=True)

    d = groups['D4']
    models = defaultdict(list)
    source_count = 0
    for w in range(d.order):
        for x in sorted(d.lower[w]):
            if d.fc[x] and d.kl(x, w) == (1, 1) and d.length[w]-d.length[x] in (3, 5):
                data = poset_data(d, x, w)
                models[(d.length[w]-d.length[x], data['fingerprint'])].append((x, w, data))
                source_count += 1
    result['search_counts']['D4_fc_bottom_P_1_plus_q_rank_3_or_5'] = source_count
    print('D4 source candidates', source_count, 'fingerprints', len(models), flush=True)
    for name in ['B3', 'B4', 'A4']:
        target = groups[name]
        counts, found = Counter(), set()
        for v in range(target.order):
            for y in sorted(target.lower[v]):
                gap = target.length[v]-target.length[y]
                if target.fc[y] or gap not in (3, 5) or target.kl(y, v) != (1, 1):
                    continue
                counts['non_fc_bottom_P_1_plus_q_rank_3_or_5'] += 1
                # A B_n support containing s_0,s_1 uses the edge labelled 4.
                # Require the full rank support, and avoid a mere constant signed suffix.
                full_support = set(target.words[v]) == set(range(target.rank))
                if not full_support:
                    continue
                counts['full_support'] += 1
                data = poset_data(target, y, v)
                options = models.get((gap, data['fingerprint']), [])
                if options:
                    counts['matching_fingerprint'] += 1
                if gap in found:
                    continue
                for x, w, source_data in options:
                    mapping = isomorphism(source_data, data)
                    if mapping is None:
                        continue
                    cert = certify(d, x, w, target, y, v, mapping)
                    cert['target_top_full_simple_support'] = full_support
                    cert['target_obstruction_to_type_A_parabolic_cosets'] = non_type_a_coset_certificate(target, y, v)
                    result['certificates'].append(cert)
                    found.add(gap)
                    print('FOUND', name, 'rank', gap, 'source', d.describe(x), d.describe(w),
                          'target', target.describe(y), target.describe(v),
                          'rank vector', cert['rank_vector'], flush=True)
                    break
        result['search_counts'][name] = dict(counts)

    # A rank-five mu=1 example uses all of the D4 fork, with P=(1+q)^2.
    mu_models = []
    for w in range(d.order):
        for x in sorted(d.lower[w]):
            if d.fc[x] and d.length[w]-d.length[x] == 5 and d.mu(x, w) == 1:
                mu_models.append((x, w, poset_data(d, x, w)))
    result['search_counts']['D4_fc_bottom_rank_5_mu_1'] = len(mu_models)
    target, target_count, found_mu = groups['B4'], 0, False
    for v in range(target.order):
        for y in sorted(target.lower[v]):
            if target.fc[y] or target.length[v]-target.length[y] != 5 or target.kl(y, v) != (1, 2, 1):
                continue
            target_count += 1
            if found_mu:
                continue
            data = poset_data(target, y, v)
            for x, w, source_data in mu_models:
                mapping = isomorphism(source_data, data)
                if mapping is None:
                    continue
                cert = certify(d, x, w, target, y, v, mapping)
                cert['target_top_full_simple_support'] = set(target.words[v]) == set(range(target.rank))
                cert['source_top_full_simple_support'] = set(d.words[w]) == set(range(d.rank))
                cert['target_obstruction_to_type_A_parabolic_cosets'] = non_type_a_coset_certificate(target, y, v)
                result['certificates'].append(cert)
                found_mu = True
                print('FOUND B4 rank5 mu1:', d.elements[x], d.elements[w], target.elements[y], target.elements[v], flush=True)
                break
    assert found_mu
    result['search_counts']['B4_non_fc_bottom_rank_5_P_1_plus_2q_plus_q2'] = target_count

    b2 = Coxeter('B', 2)
    b2_top = b2.length.index(4)
    b2_data = poset_data(b2, 0, b2_top)
    tested, matches = 0, 0
    for w in range(d.order):
        for x in d.lower[w]:
            if d.fc[x] and d.length[w]-d.length[x] == 4:
                tested += 1
                matches += int(isomorphism(b2_data, poset_data(d, x, w)) is not None)
    assert matches == 0
    result['bounded_absence_check'] = dict(
        target_type='B2', target_bottom=b2.describe(0), target_top=b2.describe(b2_top),
        target_rank_vector=b2_data['fingerprint'][0], source_type='D4',
        tested_fc_bottom_rank_four_intervals=tested, isomorphic_models=matches,
        limitation='Excludes only D4 FC-bottom models; makes no assertion about Dn for n>4.')

    certificate_path = Path(__file__).resolve().parents[1] / 'results' / 'd6-recurrence-certificate.json'
    certificate_path.parent.mkdir(exist_ok=True)
    result['independent_bad_D6_computation'] = bad_d6(certificate_path)
    result['D6_dependency_verification'] = verify(certificate_path)
    result['bad_D8_computation_status'] = 'Not attempted by the default command; the thesis base computation is inherited, not independently reproduced.'

    result['bad_element_arithmetic'] = []
    for n in range(4, 33, 2):
        ellw = Fraction(3*n*n, 8)+Fraction(n, 4)
        ellx = Fraction(n, 2)+1
        a = Fraction(3*n, 4) if n % 4 == 0 else Fraction(3*n+2, 4)
        delta = ellw-ellx
        assert all(v.denominator == 1 for v in (ellw, ellx, a, delta))
        degree = min((delta-1)//2, (ellw-a)//2)
        if n == 8:
            degree = min(degree, 9)
        result['bad_element_arithmetic'].append(dict(n=n, ell_w=int(ellw), ell_x=int(ellx),
            interval_rank=int(delta), a=int(a), degree_bound=int(degree),
            mu=1 if n == 6 else 0,
            basis='exact arithmetic in supplied thesis formulas; not a KL computation'))
    result['elapsed_seconds'] = round(perf_counter()-started, 3)
    output = Path(__file__).resolve().parents[1] / 'results' / 'computation.json'
    output.parent.mkdir(exist_ok=True)
    output.write_text(json.dumps(result, indent=2)+'\n')
    print('Saved', output, 'in', result['elapsed_seconds'], 'seconds', flush=True)


if __name__ == '__main__':
    main()
