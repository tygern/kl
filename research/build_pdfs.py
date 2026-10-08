"""Export the project's LaTeX manuscripts under results/ to PDF.

Run from any directory: python3 research/build_pdfs.py
The installed latexmk is used; no packages are installed and shell escape is
disabled. Only results/*.tex is built. The archived research notes under
research/archive/notes/ and the frozen manuscript copy inside the supplement
are not project deliverables and are deliberately skipped. Every path written
to results/pdf-build-manifest.json is repository-relative.
"""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
sources = sorted((ROOT / 'results').glob('*.tex'))
if not sources:
    raise SystemExit('No LaTeX sources found under results/.')
if len({p.stem for p in sources}) != len(sources):
    raise SystemExit('Duplicate TeX basenames need distinct output names.')
compiler = shutil.which('latexmk')
if not compiler:
    raise SystemExit('The installed latexmk executable is unavailable.')
out = ROOT / 'output' / 'pdf'
out.mkdir(parents=True, exist_ok=True)


def rel(path):
    """Repository-relative POSIX path for manifests."""
    return Path(path).resolve().relative_to(ROOT).as_posix()


records = []
log = []
for source in sources:
    build = ROOT / 'tmp' / 'pdfs' / 'build' / source.stem
    build.mkdir(parents=True, exist_ok=True)
    # Run from the repository root with relative paths, so the recorded
    # command is exactly the command that ran and contains no machine paths.
    command = ['latexmk', '-pdf', '-interaction=nonstopmode', '-halt-on-error',
               '-file-line-error', '-no-shell-escape', f'-outdir={rel(build)}', rel(source)]
    result = subprocess.run([compiler] + command[1:], cwd=ROOT, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    log.append(result.stdout)
    (ROOT / 'results' / 'pdf-build.log').write_text('\n'.join(log))
    if result.returncode:
        raise SystemExit(f'Build failed: {rel(source)}; see results/pdf-build.log')
    generated = build / f'{source.stem}.pdf'
    data = generated.read_bytes()
    assert data.startswith(b'%PDF-')
    target = out / generated.name
    shutil.copy2(generated, target)
    records.append(dict(source=rel(source), pdf=rel(target),
                        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                        pdf_sha256=hashlib.sha256(data).hexdigest(), bytes=len(data),
                        command=command, cwd='repository root'))
(ROOT / 'results' / 'pdf-build-manifest.json').write_text(json.dumps(records, indent=2) + '\n')
for row in records:
    print(row['pdf'])
