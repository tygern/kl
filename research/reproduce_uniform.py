"""Reproduce rank-uniform certificates with Python's standard library only.

The mathematical all-rank proof is in en_uniform/construction.txt and
en_uniform/referee.txt. Arithmetic samples alone are not that proof.
Run: uv run --offline --no-project --python /opt/homebrew/bin/python3 research/reproduce_uniform.py
"""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT / 'research' / 'en_uniform'
log = []
for name in ('construction.py', 'referee_verify.py'):
    result = subprocess.run([sys.executable, str(HERE/name)], cwd=ROOT,
                            text=True, capture_output=True)
    log.append(f'{name}\n{result.stdout}{result.stderr}')
    (ROOT/'results'/'uniform-reproduction.log').write_text('\n'.join(log))
    if result.returncode:
        raise SystemExit(f'{name} failed; inspect results/uniform-reproduction.log')
    print(f'{name}: passed', flush=True)
files = [p for p in HERE.iterdir() if p.suffix in {'.py', '.json', '.txt'}]
files.extend((Path(__file__).resolve(), ROOT/'results'/'exceptional-leading.tex'))
manifest = dict(status='both exact certificate runners passed',
                scope='Symbolic identities plus finite arithmetic counterchecks; all-r proof separately audited.',
                files={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in sorted(files)})
(ROOT/'results'/'uniform-reproduction-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
