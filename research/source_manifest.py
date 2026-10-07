"""Record downloaded public sources with immutable local SHA-256 digests."""
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
URLS = {
    'combinatorial-invariance-2026.pdf': 'https://raw.githubusercontent.com/openai/math/main/preprints/Combinatorial-Invariance-of-Kazhdan-Lusztig-Polynomials-September-24-2026/paper.pdf',
    'gern-thesis-2013.pdf': 'https://arxiv.org/pdf/1304.6074',
    'green-leading-2008.pdf': 'https://arxiv.org/pdf/0801.1650',
    'green-jones-traces-2007.pdf': 'https://arxiv.org/pdf/math/0509362',
    'green-losonczy-cells-2001.pdf': 'https://arxiv.org/pdf/math/0102003',
    'jones-deodhar-2007.pdf': 'https://arxiv.org/pdf/0711.1391',
    'woo-patterns-2006.pdf': 'https://arxiv.org/pdf/math/0611328',
}
records = []
for filename, url in URLS.items():
    path = ROOT / 'sources' / filename
    data = path.read_bytes()
    if not data.startswith(b'%PDF-'):
        raise ValueError(f'Not a PDF: {path}')
    records.append(dict(file=f'sources/{filename}', source_url=url,
                        size_bytes=len(data), sha256=sha256(data).hexdigest()))
manifest = dict(research_date_local='2026-10-06', timezone='America/Chicago',
                note='Locally archived snapshots, not a guarantee of external acceptance or later version stability.',
                sources=records,
                not_archived=[dict(title='Chmutov 2014 thesis',
                    url='https://hdl.handle.net/2027.42/108802',
                    reason='University PDF returned an access challenge; indexed primary text and author-uploaded text inspected; see literature audit.')])
(ROOT / 'sources' / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
print(f'Archived source manifest: {len(records)} verified PDF files.')
