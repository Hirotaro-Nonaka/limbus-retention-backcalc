"""Render the paper drafts (markdown) to submission-ready PDFs.

Pipeline: markdown -> styled HTML (print CSS, CJK-aware fonts) -> PDF via
Microsoft Edge headless printing. Run from the repo root:

    py scripts/make_pdf.py

Outputs paper/Limbus_Retention_Paper_JA.pdf and paper/Limbus_Retention_Paper_EN.pdf.
"""
import os, subprocess, sys, tempfile, time
from pathlib import Path

import markdown

ROOT = Path(__file__).resolve().parent.parent
PAPER = ROOT / 'paper'
CHROME = r'C:\Program Files\Google\Chrome\Application\chrome.exe'
EDGE = r'C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe'
# Prefer Chrome: on this machine Edge 150 headless intermittently exits 0
# without printing (see repo history); Chrome prints reliably.
BROWSER = CHROME if Path(CHROME).exists() else EDGE

CSS = """
@page { size: A4; margin: 22mm 20mm; }
body {
  font-family: {FONTS};
  font-size: 10.5pt; line-height: 1.7; color: #111;
  max-width: 170mm; margin: 0 auto;
}
h1 { font-size: 16pt; line-height: 1.4; border-bottom: 2px solid #333; padding-bottom: 6px; }
h2 { font-size: 13pt; margin-top: 1.6em; border-bottom: 1px solid #999; padding-bottom: 3px; }
h3 { font-size: 11.5pt; margin-top: 1.3em; }
h4 { font-size: 10.5pt; margin-top: 1.2em; }
table { border-collapse: collapse; margin: 0.8em 0; font-size: 9pt; width: 100%; }
th, td { border: 1px solid #888; padding: 3px 6px; text-align: left; vertical-align: top; }
th { background: #f0f0f0; }
img { max-width: 100%; height: auto; display: block; margin: 0.8em auto; }
strong { font-weight: bold; }
blockquote { border-left: 3px solid #bbb; margin-left: 0; padding-left: 1em; color: #333; }
hr { border: none; border-top: 1px solid #bbb; margin: 1.5em 0; }
li { margin: 0.25em 0; }
h1, h2, h3, h4 { page-break-after: avoid; }
table, img { page-break-inside: avoid; }
"""

JP_FONTS = '"Yu Mincho", "Hiragino Mincho ProN", "MS Mincho", serif'
EN_FONTS = 'Georgia, "Times New Roman", serif'


def build(md_name: str, pdf_name: str, fonts: str, lang: str):
    src = PAPER / md_name
    text = src.read_text(encoding='utf-8')
    body = markdown.markdown(text, extensions=['tables', 'fenced_code'])
    html = (
        f'<!DOCTYPE html><html lang="{lang}"><head><meta charset="utf-8">'
        f'<style>{CSS.replace("{FONTS}", fonts)}</style></head>'
        f'<body>{body}</body></html>'
    )
    html_path = PAPER / (Path(pdf_name).stem + '.html')
    html_path.write_text(html, encoding='utf-8')
    pdf_path = PAPER / pdf_name
    # Isolated profile + --no-first-run: without these, a resident Edge/WebView2
    # process can take over the launch and silently skip printing (exit code 0,
    # no output file).
    profile = Path(tempfile.gettempdir()) / 'make_pdf_edge_profile'
    tmp_pdf = pdf_path.with_name('_build_' + pdf_path.name)
    tmp_pdf.unlink(missing_ok=True)
    cmd = [BROWSER, f'--user-data-dir={profile}', '--no-first-run', '--no-default-browser-check',
           '--headless=new', '--disable-gpu', '--no-pdf-header-footer',
           f'--print-to-pdf={tmp_pdf}', html_path.resolve().as_uri()]
    # Edge's first headless launch sometimes exits 0 without printing; retry.
    for attempt in range(4):
        subprocess.run(cmd, check=True, timeout=120)
        if tmp_pdf.exists():
            break
        time.sleep(2)
    else:
        raise RuntimeError(f'browser never wrote {tmp_pdf} after retries')
    os.replace(tmp_pdf, pdf_path)
    print(f'wrote {pdf_path}')


if __name__ == '__main__':
    build('Limbus_Retention_Paper_draft.md', 'Limbus_Retention_Paper_JA.pdf', JP_FONTS, 'ja')
    build('Limbus_Retention_Paper_draft_EN.md', 'Limbus_Retention_Paper_EN.pdf', EN_FONTS, 'en')
