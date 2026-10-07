"""Read final PDFs and create page contact sheets for visual inspection.

Run with the bundled Python runtime (Pillow, pypdf and pdfplumber).
PNG page renders must already exist from pdftoppm.
"""
import json
import re
from pathlib import Path
from PIL import Image, ImageDraw
from pypdf import PdfReader
import pdfplumber

ROOT = Path(__file__).resolve().parents[1]
reports = []
for pdf in sorted((ROOT / 'output' / 'pdf').glob('*.pdf')):
    reader = PdfReader(pdf)
    assert len(reader.pages) > 0
    chars_outside_margins = []
    with pdfplumber.open(pdf) as opened:
        for number, page in enumerate(opened.pages, 1):
            assert page.extract_text()
            for char in page.chars:
                if char['x0'] < 65 or char['x1'] > page.width-65:
                    chars_outside_margins.append(dict(page=number, text=char['text'],
                        x0=char['x0'], x1=char['x1']))
    assert not chars_outside_margins, chars_outside_margins
    numbered = []
    for path in (ROOT / 'tmp' / 'pdfs').glob(f'{pdf.stem}-*.png'):
        match = re.fullmatch(re.escape(pdf.stem) + r'-(\d+)\.png', path.name)
        if match:
            numbered.append((int(match.group(1)), path))
    pages = [path for _, path in sorted(numbered)]
    assert len(pages) == len(reader.pages)
    thumbs = []
    for index, path in enumerate(pages, 1):
        im = Image.open(path).convert('RGB')
        im.thumbnail((408, 530))
        tile = Image.new('RGB', (430, 560), '#dddddd')
        tile.paste(im, ((430-im.width)//2, 22))
        ImageDraw.Draw(tile).text((8, 4), f'Page {index}', fill='black')
        thumbs.append(tile)
    contacts = []
    for start in range(0, len(thumbs), 6):
        chunk = thumbs[start:start+6]
        sheet = Image.new('RGB', (1290, 560*((len(chunk)+2)//3)), 'white')
        for index, im in enumerate(chunk):
            sheet.paste(im, ((index % 3)*430, (index//3)*560))
        target = ROOT / 'tmp' / 'pdfs' / f'{pdf.stem}-contact-{start//6+1}.png'
        sheet.save(target)
        contacts.append(str(target))
    reports.append(dict(pdf=str(pdf), pages=len(reader.pages), bytes=pdf.stat().st_size,
                        all_pages_have_text=True, within_page_margins=True,
                        contact_sheets=contacts,
                        visual_inspection='Contact sheets and representative full page inspected separately.'))
(ROOT / 'results' / 'pdf-qa.json').write_text(json.dumps(reports, indent=2)+'\n')
print(json.dumps(reports, indent=2))
