#!/usr/bin/env python3
"""Reescala um PDF de uma página para a largura final (mm), como vetor: o conteúdo e a página são multiplicados pelo mesmo
fator; fontes embutidas são mantidas. Uso: scale_pdf.py entrada.pdf saida.pdf largura_mm  (usado por render_final_pdf.R)"""
import sys
from pypdf import PdfReader, PdfWriter

src, dst, width_mm = sys.argv[1], sys.argv[2], float(sys.argv[3])
r = PdfReader(src)
assert len(r.pages) == 1
page = r.pages[0]
w_pt = float(page.mediabox.width)
f = (width_mm / 25.4 * 72) / w_pt
page.scale_by(f)
w = PdfWriter(); w.add_page(page)
with open(dst, 'wb') as fh: w.write(fh)
print(f'escala {f:.4f}; página final {float(page.mediabox.width) / 72 * 25.4:.1f} x {float(page.mediabox.height) / 72 * 25.4:.1f} mm')
