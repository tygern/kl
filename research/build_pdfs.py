"""Export every project LaTeX source to PDF using the installed TeX toolchain.

The native editor compiler is used separately for document diagnostics.
Run from any directory: python3 research/build_pdfs.py
No packages are installed and shell escape is disabled.
"""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
# The supplement contains a frozen copy of the manuscript, not a fourth
# editable document. Exclude archive staging just as we exclude PDF builds.
sources = sorted(p for p in ROOT.rglob('*.tex')
                 if not any(part in {'.git', 'tmp', 'output', 'supplement'} for part in p.relative_to(ROOT).parts))
if not sources:
    raise SystemExit('No project LaTeX sources found.')
if len({p.stem for p in sources}) != len(sources):
    raise SystemExit('Duplicate TeX basenames need distinct output names.')
compiler = shutil.which('latexmk')
if not compiler:
    raise SystemExit('The installed latexmk executable is unavailable.')
out = ROOT / 'output' / 'pdf'
out.mkdir(parents=True, exist_ok=True)
records = []
log = []
for source in sources:
    build = ROOT / 'tmp' / 'pdfs' / 'build' / source.stem
    build.mkdir(parents=True, exist_ok=True)
    command = [compiler, '-pdf', '-interaction=nonstopmode', '-halt-on-error',
               '-file-line-error', '-no-shell-escape', f'-outdir={build}', str(source)]
    result = subprocess.run(command, cwd=source.parent, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    log.append(result.stdout)
    (ROOT / 'results' / 'pdf-build.log').write_text('\n'.join(log))
    if result.returncode:
        raise SystemExit(f'Build failed: {source}; see results/pdf-build.log')
    generated = build / f'{source.stem}.pdf'
    data = generated.read_bytes()
    assert data.startswith(b'%PDF-')
    target = out / generated.name
    shutil.copy2(generated, target)
    records.append(dict(source=str(source), pdf=str(target),
                        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                        pdf_sha256=hashlib.sha256(data).hexdigest(), bytes=len(data),
                        command=command))
(ROOT / 'results' / 'pdf-build-manifest.json').write_text(json.dumps(records, indent=2)+'\n')
for row in records:
    print(row['pdf'])
