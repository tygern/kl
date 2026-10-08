#!/usr/bin/env python3
"""Rebuild and verify the local proof supplement; Python standard library only.

Default mode recomputes every classification, catalogue and certificate that
the manuscript's proofs depend on, verifies the shipped certificates against
the regenerated ones, and compares the stable mathematical outcomes with
expected-summary.json.  --full additionally reruns the historical full-group
E6/E7 searches, the second (geometric) model of the D6/D8 Kazhdan-Lusztig
computation and the larger uniform-family checks, comparing every
regenerated snapshot with the shipped one.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parent
# Generated JSON files from the finite terminal engines carry wall-clock
# timing; these keys are ignored wherever snapshots are compared.
VOLATILE_KEYS = {'seconds', 'elapsed_seconds'}


def strip_volatile(value):
    if isinstance(value, dict):
        return {k: strip_volatile(v) for k, v in value.items() if k not in VOLATILE_KEYS}
    if isinstance(value, list):
        return [strip_volatile(v) for v in value]
    return value


def canonical(path):
    return strip_volatile(json.loads(Path(path).read_text()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full',action='store_true',help='also rerun the historical E6/E7 full-group searches, the geometric D6/D8 model and the larger uniform-family checks')
    parser.add_argument('--compiler',help='installed C++17 compiler executable; default searches clang++, c++, g++')
    args = parser.parse_args()
    if sys.version_info < (3,10): raise SystemExit('Python 3.10 or later is required.')
    if not __debug__: raise SystemExit('Do not use python -O: proof assertions must remain enabled.')
    manifest = json.loads((ROOT/'MANIFEST.json').read_text())
    for name,digest in manifest['sha256'].items():
        file = ROOT/name
        if not file.is_file() or hashlib.sha256(file.read_bytes()).hexdigest()!=digest:
            raise SystemExit('Input checksum mismatch: '+name)
    inputs = json.loads((ROOT/'INPUTS.json').read_text())
    manuscript_sha256_at_build = next(f['source_sha256'] for f in inputs['files'] if f['source_path']=='results/exceptional-leading.tex')
    compiler = args.compiler or next((shutil.which(c) for c in ('clang++','c++','g++') if shutil.which(c)),None)
    if not compiler: raise SystemExit('An installed C++17 compiler is required; nothing will be installed.')
    (ROOT/'runs').mkdir(exist_ok=True)
    run = Path(tempfile.mkdtemp(prefix='full-' if args.full else 'default-',dir=ROOT/'runs'))
    work = run/'work'
    shutil.copytree(ROOT/'payload',work)
    logs=run/'logs';logs.mkdir()
    bins=run/'bin';bins.mkdir()
    started=time.monotonic()
    steps=[]
    snapshots=[]
    env=os.environ.copy();env.pop('PYTHONOPTIMIZE',None)
    env['PYTHONDONTWRITEBYTECODE']='1'
    # Avoid inherited Python import paths affecting the supplement's modules.
    env.pop('PYTHONPATH',None)
    def execute(name,command,output=None):
        before=time.monotonic()
        with (logs/(name+'.stdout.log')).open('w') as out, (logs/(name+'.stderr.log')).open('w') as err:
            p=subprocess.run(list(map(str,command)),cwd=work,env=env,stdout=out,stderr=err)
        steps.append({'step':name,'command':[str(c).replace(str(run),'<run>') for c in command], 'returncode':p.returncode,
                      'seconds':round(time.monotonic()-before,3)})
        (run/'steps.json').write_text(json.dumps(steps,indent=2)+'\n')
        if p.returncode: raise RuntimeError(f'{name} failed; inspect {logs}')
        if output:
            target=work/output;target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(logs/(name+'.stdout.log'),target)
        print(f'{name}: passed',flush=True)
    def build(name,source):
        exe=bins/name
        execute('build-'+name,[compiler,'-O3','-std=c++17',source,'-o',exe])
        return exe
    # -I omits the script directory, so local-import scripts require the ordinary
    # interpreter with sanitized PYTHONPATH. These trusted inputs are hashed.
    def py(name,path,*extra):execute(name,[sys.executable,path,*extra])
    def read(path):return json.loads((work/path).read_text())
    def compare_snapshot(path,select=lambda d:d):
        """The regenerated file must equal the shipped one up to volatile timing keys."""
        regenerated=select(canonical(work/path)); shipped=select(canonical(ROOT/'payload'/path))
        if regenerated!=shipped: raise AssertionError('Regenerated '+path+' differs from the shipped snapshot; inspect '+str(work/path))
        snapshots.append(path); print(f'snapshot {path}: regenerated output equals the shipped file',flush=True)
    try:
        flat=build('flat','research/en_e8/parabolic_terminals.cpp')
        recursive=build('recursive','research/en_e8/recursive_terminals.cpp')
        d7=build('d7','research/en_independent/e8_d7_cosets.cpp')
        fc=build('fc','research/en_independent/fc_catalogue.cpp')
        e6mu=build('e6-mu-table','research/review_checks/e6_mu_table.cpp')
        e9kl=build('e9-quotient-kl','research/review_checks/e9_quotient_kl.cpp')
        d8kl=build('d8-gern-kl','research/review_checks/d8_gern_kl.cpp')
        fcmax=build('fc-maxima','research/review_checks/fc_maxima.cpp')
        for n in (6,7,8,9):
            execute(f'E{n}-FC',[fc,n],f'research/en_independent/e{n}-fc.json')
        for n in (6,7,8):
            flatfile=f'research/en_e8/e{n}-validation.json' if n<8 else 'research/en_e8/e8-terminals.json'
            execute(f'E{n}-flat',[flat,n],flatfile)
            execute(f'E{n}-recursive',[recursive,n],f'research/en_e8/e{n}-recursive.json')
        execute('E8-D7',[d7],'research/en_independent/e8-d7-terminals.json')
        if args.full:
            rootids=build('historical-rootids','research/broad_exceptional/enumerate_bad.cpp')
            matrix=build('historical-e7-matrix','research/verify_e7_matrices.cpp')
            for n in (6,7):
                execute(f'E{n}-historical-rootids',[rootids,n],f'research/broad_exceptional/e{n}_bad.json')
                compare_snapshot(f'research/broad_exceptional/e{n}_bad.json')
            py('E6-historical-matrix','research/verify_e6_independent.py')
            compare_snapshot('results/e6-independent-certificate.json')
            execute('E7-historical-matrix',[matrix],'results/e7-independent-certificate.json')
            compare_snapshot('results/e7-independent-certificate.json')
        py('finite-source-snapshot-audit','research/verify_exceptional.py')
        py('finite-method-comparison','research/en_e8/verify_outputs.py')
        py('E8-independent-verification','research/en_independent/verify_e8.py')
        py('finite-support-matchings','research/ai-review-notes/check_finite_descents.py')
        compare_snapshot('research/ai-review-notes/finite-descents.json')
        py('D6-recurrence-certificate','computations/verify_certificate.py')
        py('uniform-symbolic-and-matrix-checks','research/en_uniform/referee_verify.py')
        # The uniform verifier records the hash of the manuscript it audited; it
        # must be the manuscript snapshot recorded by the builder in INPUTS.json.
        audited=read('research/en_uniform/referee-certificate.json')['audited_source_sha256']
        if audited!=manuscript_sha256_at_build: raise AssertionError('Uniform verifier audited a manuscript other than the one recorded at build time.')
        py('affine-proved-family','affine_proof.py')
        compare_snapshot('research/en_affine_referee/proved-affine-certificate.json')
        py('affine-FC-closure-and-covers','research/en_affine_referee/verify_fc_catalogue.py')
        py('E10-all-k','research/en_affine_referee/verify_indefinite_e10.py')
        # Certificates added after the review of 7 October 2026 (Section 5 of the manuscript
        # and the remarks of Sections 3 and 4). Each program asserts its expected values
        # itself; the runner additionally requires equality with the shipped certificate.
        execute('E6-complete-mu-table',[e6mu,'E','6','results/e6-mu-table.json'])
        compare_snapshot('results/e6-mu-table.json')
        execute('E9-odd-gap-polynomials',[e9kl,'certify','results/e9-odd-gap-certificate.json'])
        compare_snapshot('results/e9-odd-gap-certificate.json')
        execute('D6-D8-Gern-polynomials',[d8kl,'results/d8-gern-certificate.json','0','both' if args.full else 'sp'])
        compare_snapshot('results/d8-gern-certificate.json',lambda d:(d['ranks'],d['status']))
        execute('FC-maxima-E6-E13',[fcmax,'results/fc-maxima-certificate.json',6,7,8,9,10,11,12,13])
        compare_snapshot('results/fc-maxima-certificate.json')
        py('terminal-structure','research/review_checks/terminal_structure.py')
        compare_snapshot('results/terminal-structure-certificate.json')
        py('uniform-family-covers','research/review_checks/uniform_family_checks.py',*(['--full'] if args.full else []))
        def uniform_rows(d):
            return ({(r['r'],r['k']):r for r in d['pairs']},{(c['r'],c['k']):c for c in d['bruhat_covers']},d['status'])
        if args.full:
            # The full run checks more pairs; the shipped (default) rows must reappear unchanged.
            regenerated=uniform_rows(canonical(work/'results/uniform-family-certificate.json'))
            shipped=uniform_rows(canonical(ROOT/'payload/results/uniform-family-certificate.json'))
            if not(regenerated[2]==shipped[2]=='passed' and all(regenerated[i][k]==v for i in (0,1) for k,v in shipped[i].items())):
                raise AssertionError('Full uniform-family run disagrees with the shipped certificate.')
            snapshots.append('results/uniform-family-certificate.json (shipped rows)')
        else:
            compare_snapshot('results/uniform-family-certificate.json')
        py('affine-D4-mu-2','research/verify_affine_d4_r.py')
        compare_snapshot('results/affine-d4-independent.json')
        summary={'finite':{}}
        for n in (6,7,8):
            row=read(f'research/en_e8/e{n}-recursive.json')
            summary['finite'][f'E{n}']={k:row[k] for k in ('group_order','right_terminal_count','commuting_terminals','noncommuting_terminal_count','fc_count')}
            summary['finite'][f'E{n}']['bad_lengths']=sorted(x['length'] for x in row['bad'])
        d=read('research/en_independent/e8-d7-terminals.json')
        summary['E8_D7']={k:d[k] for k in ('cosets','parabolic_right_terminals','candidates_tested','terminal_count')}
        d=json.loads((logs/'D6-recurrence-certificate.stdout.log').read_text())
        summary['D6']={k:d[k] for k in ('records','root_polynomial','root_right_descent_checks','root_R_reciprocity','evaluator_called')}
        d=read('research/en_uniform/referee-certificate.json')
        summary['uniform']={'symbolic_assertions':d['symbolic']['all_assertions_pass'],'checked_ranks':[r['r'] for r in d['rows']]}
        d=read('research/en_affine_referee/proved-affine-certificate.json')
        summary['affine']={'length':[d['length_formula']['intercept'],d['length_formula']['slope']], 'terminal_for_all_k':d['terminal_for_all_k']}
        d=read('research/en_affine_referee/fc-cover-certificate.json')
        summary['affine'].update({k:d[k] for k in ('FC_count','maximum_FC_length','base_Bruhat_cover_count','FC_base_covers')})
        d=read('research/en_affine_referee/e10-certificate.json')
        summary['E10']={k:d[k] for k in ('rank','base_length','T_minus_identity_cube_zero','full_support_for_all_k','maximum_independent_size','eligible_mu')}
        matching=read('research/ai-review-notes/finite-descents.json')
        summary['finite_matchings']=matching
        d=read('results/e6-mu-table.json')
        summary['E6_mu']={k:d[k] for k in ('group_order','fully_commutative_count','max_mu','mu_histogram','pairs_with_mu_at_least_2',
                                            'fully_commutative_lower_endpoint','fully_commutative_upper_endpoint','pairs_attaining_max_mu',
                                            'largest_mu_with_fully_commutative_upper_endpoint')}
        d=read('results/e9-odd-gap-certificate.json')
        summary['E9_odd_gap']=[{k:e[k] for k in ('name','word_string','length','L','R','involution','terminal','quotient_size','lower_ideal_size','eligible_fully_commutative_bottoms')} for e in d['elements']]
        d=read('results/d8-gern-certificate.json')
        summary['D6_D8_Gern']={rank:{k:v[k] for k in ('w_length','x_length','lower_ideal_size','interval_size','P_x_w_ascending','P_e_w_ascending','mu_x_w')} for rank,v in d['ranks'].items()}
        d=read('results/fc-maxima-certificate.json')
        summary['FC_maxima']={rank:{k:v[k] for k in ('fully_commutative_count','maximum_length')} for rank,v in d['ranks'].items()}
        d=read('results/terminal-structure-certificate.json')
        summary['terminal_structure']={'rows_not_from_type_D':{k:v['rows_not_from_type_D'] for k,v in d['gern_plus_two'].items() if isinstance(v,dict)},
                                       'type_D_noncommuting_terminals':{k:v['lengths'] for k,v in d['type_D_enumeration'].items()},
                                       'negated_roots':{k:v['negated_positive_roots'] for k,v in d['structure'].items()},
                                       'layers':{k:v['layers'] for k,v in d['structure'].items()},
                                       'chain_E8':[[r['lower'],r['upper'],r['prefix'] and r['suffix']] for r in d['right_weak_order_chain']['E8']['relations']],
                                       'inner_longest_parabolic_lengths':{k:v['max_inner_longest_parabolic_length'] for k,v in d['inner_longest_parabolic_factors'].items()},
                                       'w0_factorizations':{k:v['w0_factorization'] for k,v in d['inner_longest_parabolic_factors'].items() if 'w0_factorization' in v}}
        d=canonical(ROOT/'payload/results/uniform-family-certificate.json')  # default rows, identical in both modes (checked above)
        summary['uniform_family']={'checked_pairs':d['checked_pairs'],'lengths':[r['length'] for r in d['pairs']],
                                   'length_formula_holds':d['length_formula_holds_for_all_checked_pairs'],
                                   'bruhat_covers':[[c['r'],c['k'],c['distinct_covers'],c['fully_commutative_covers']] for c in d['bruhat_covers']]}
        d=read('results/affine-d4-independent.json')
        summary['affine_D4']={k:d[k] for k in ('lengths','lower_ideal_size','polynomial','mu')}
        expected=json.loads((ROOT/'expected-summary.json').read_text())
        (run/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
        if summary != expected:raise AssertionError('Stable proof outputs differ from expected-summary.json; inspect '+str(run/'summary.json'))
        # Also match the table words checked by the independent matching verifier
        # to the actual generated terminal matrices, not merely their lengths.
        import importlib.util
        spec=importlib.util.spec_from_file_location('comparison',work/'research/en_e8/verify_outputs.py')
        comparison=importlib.util.module_from_spec(spec);spec.loader.exec_module(comparison)
        # check_certificate uses module ROOT, which is the isolated working tree.
        for n in (6,7,8):
            _, generated=comparison.check_certificate(f'research/en_e8/e{n}-recursive.json',n)
            table=[row for row in matching if row['type']==f'E{n}']
            edges=[(j,j+1) for j in range(n-2)]+[(2,n-1)]
            def matrix(word):
                cols=tuple(tuple(int(i==j) for i in range(n)) for j in range(n))
                for s in map(int,word):
                    old=cols;out=list(old);out[s]=tuple(-v for v in old[s])
                    for a,b in edges:
                        if s in (a,b):
                            t=b if s==a else a
                            out[t]=tuple(u+v for u,v in zip(old[t],old[s]))
                    cols=tuple(out)
                return cols
            assert {matrix(row['word']) for row in table} == set(generated)
        status={'status':'passed','mode':'full' if args.full else 'default',
                'seconds':round(time.monotonic()-started,3),'compiler':str(compiler),
                'python':sys.version,'steps':len(steps),'work':str(work.relative_to(ROOT)),
                'expected_outputs_match':True,'printed_table_matches_generated_terminals':True,
                'manuscript_sha256_matches_build_record':True,'snapshots_compared':snapshots}
        (run/'status.json').write_text(json.dumps(status,indent=2)+'\n')
        print('All proof checks and expected outputs passed. Logs and certificates: '+str(run),flush=True)
    except Exception:
        print('Verification failed. Retained working tree and logs: '+str(run),file=sys.stderr)
        raise

if __name__=='__main__':main()
