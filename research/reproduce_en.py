"""Rebuild the proved finite/infinite E_n certificates; standard library only.

Requires an installed C++17 compiler. No network or package installation.
Exploratory random searches and capped lower-ideal probes are not rerun.
Run: python3 research/reproduce_en.py
"""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
LOGS = ROOT / 'research/en_independent/reproduction-logs'

def run(command, name, output=None):
    result = subprocess.run(command, cwd=ROOT, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (LOGS / (name + '.stdout.log')).write_text(result.stdout)
    (LOGS / (name + '.stderr.log')).write_text(result.stderr)
    if result.returncode:
        raise RuntimeError(f'{name} failed; inspect {LOGS}')
    if output:
        (ROOT / output).write_text(result.stdout)
    print(name + ': passed', flush=True)

def main():
    LOGS.mkdir(parents=True, exist_ok=True)
    compiler = shutil.which('clang++') or shutil.which('c++')
    if not compiler:
        raise SystemExit('An installed C++17 compiler is required.')
    with tempfile.TemporaryDirectory(prefix='kl-en-certificates-') as tmp:
        fc = str(Path(tmp) / 'fc')
        d7 = str(Path(tmp) / 'd7')
        run([compiler, '-O3', '-std=c++17', 'research/en_independent/fc_catalogue.cpp', '-o', fc], 'build-fc')
        run([compiler, '-O3', '-std=c++17', 'research/en_independent/e8_d7_cosets.cpp', '-o', d7], 'build-d7')
        for n in range(6, 11):
            run([fc, str(n)], f'E{n}-FC', f'research/en_independent/e{n}-fc.json')
        run([d7], 'E8-D7-terminal', 'research/en_independent/e8-d7-terminals.json')
        for path in [
            'research/en_e8/reproduce.py',
            'research/en_independent/verify_e8.py',
            'research/en_families/positive_affine_family.py',
            'research/en_families/indefinite_positive_families.py',
            'research/en_affine_referee/verify_families.py',
            'research/en_affine_referee/verify_fc_catalogue.py',
            'research/en_affine_referee/verify_indefinite.py',
            'research/en_affine_referee/verify_indefinite_e10.py',
        ]:
            run([sys.executable, path], Path(path).parent.name + '-' + Path(path).stem)
    files = []
    for folder in ['research/en_independent', 'research/en_e8', 'research/en_families', 'research/en_affine_referee']:
        files.extend(p for p in (ROOT / folder).iterdir() if p.suffix in {'.py','.cpp','.json'})
    files.append(ROOT / 'results/e8-independent-audit.json')
    manifest = {'status': 'all proof certificate reproductions passed',
                'excluded': 'exploratory random searches and capped unresolved KL probes',
                'files': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)}}
    (ROOT / 'results/en-reproduction-manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')

if __name__ == '__main__':
    main()
