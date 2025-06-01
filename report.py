import pandas as pd, io
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter

def to_csv(muts: list[dict]) -> bytes:
    return pd.DataFrame(muts).to_csv(index=False).encode()

def to_pdf(muts: list[dict]) -> io.BytesIO:
    buf = io.BytesIO(); c = canvas.Canvas(buf, pagesize=letter)
    txt = c.beginText(40, 750)
    txt.setFont("Helvetica", 12)
    if not muts:
        txt.textLine("Tidak ada mutasi.")
    else:
        for m in muts:
            txt.textLine(f"Pos {m['Position']}: {m['Mutation']} – {m['DrugClass']} ({m['Interpretation']})")
    c.drawText(txt); c.showPage(); c.save(); buf.seek(0)
    return buf
