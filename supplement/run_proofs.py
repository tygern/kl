#!/usr/bin/env python3
"""Rebuild and verify the local proof supplement; Python standard library only."""
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full',action='store_true',help='also rebuild historical E6/E7 full-group root-ID and integer-matrix searches')
    parser.add_argument('--compiler',help='installed C++17 compiler executable; default searches clang++, c++, g++')
    args = parser.parse_args()
    if sys.version_info < (3,10): raise SystemExit('Python 3.10 or later is required.')
    if not __debug__: raise SystemExit('Do not use python -O: proof assertions must remain enabled.')
    manifest = json.loads((ROOT/'MANIFEST.json').read_text())
    for name,digest in manifest['sha256'].items():
        file = ROOT/name
        if not file.is_file() or hashlib.sha256(file.read_bytes()).hexdigest()!=digest:
            raise SystemExit('Input checksum mismatch: '+name)
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
    env=os.environ.copy();env.pop('PYTHONOPTIMIZE',None)
    env['PYTHONDONTWRITEBYTECODE']='1'
    # Avoid inherited Python import paths affecting the supplement's modules.
    env.pop('PYTHONPATH',None)
    def execute(name,command,output=None):
        before=time.monotonic()
        with (logs/(name+'.stdout.log')).open('w') as out, (logs/(name+'.stderr.log')).open('w') as err:
            p=subprocess.run(list(map(str,command)),cwd=work,env=env,stdout=out,stderr=err)
        steps.append({'step':name,'command':list(map(str,command)), 'returncode':p.returncode,
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
    def py(name,path):execute(name,[sys.executable,path])
    try:
        flat=build('flat','research/en_e8/parabolic_terminals.cpp')
        recursive=build('recursive','research/en_e8/recursive_terminals.cpp')
        d7=build('d7','research/en_independent/e8_d7_cosets.cpp')
        fc=build('fc','research/en_independent/fc_catalogue.cpp')
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
            py('E6-historical-matrix','research/verify_e6_independent.py')
            execute('E7-historical-matrix',[matrix],'results/e7-independent-certificate.json')
        py('finite-source-snapshot-audit','research/verify_exceptional.py')
        py('finite-method-comparison','research/en_e8/verify_outputs.py')
        py('E8-independent-verification','research/en_independent/verify_e8.py')
        py('finite-support-matchings','research/journal-review/check_finite_descents.py')
        py('D6-recurrence-certificate','computations/verify_certificate.py')
        py('uniform-symbolic-and-matrix-checks','research/en_uniform/referee_verify.py')
        py('affine-proved-family','affine_proof.py')
        py('affine-FC-closure-and-covers','research/en_affine_referee/verify_fc_catalogue.py')
        py('E10-all-k','research/en_affine_referee/verify_indefinite_e10.py')
        def read(path):return json.loads((work/path).read_text())
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
        matching=read('research/journal-review/finite-descents.json')
        summary['finite_matchings']=matching
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
            # Existing independent comparison verifies the supplied table-based
            # E6/E7 snapshots; E8 table words are checked explicitly below.
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
                'python':sys.version,'steps':len(steps),'work':str(work),
                'expected_outputs_match':True,'printed_table_matches_generated_terminals':True}
        (run/'status.json').write_text(json.dumps(status,indent=2)+'\n')
        print('All proof checks and expected outputs passed. Logs and certificates: '+str(run),flush=True)
    except Exception:
        print('Verification failed. Retained working tree and logs: '+str(run),file=sys.stderr)
        raise

if __name__=='__main__':main()
