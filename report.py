import pandas as pd, io
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter

def to_csv(muts: list[dict]) -> bytes:
    return pd.DataFrame(muts).to_csv(index=False).encode()

def to_pdf(muts: list[dict]) -> io.BytesIO:
    """
    Buat PDF dengan format per‐entry seperti:

    Position    : ...
    Mutation    : ...
    Drug Class  : ...
    Interpretasi:
        <poin pertama>
        <poin kedua>
        ...
    Sample ID   : ...
    Reference   : ...

    (Kosongkan satu baris, lalu entry berikutnya)
    """
    buf = io.BytesIO()
    c = canvas.Canvas(buf, pagesize=letter)
    text_object = c.beginText(40, 750)
    text_object.setFont("Helvetica", 11)

    if not muts:
        text_object.textLine("Tidak ada mutasi.")
    else:
        for idx, m in enumerate(muts):
            posisi = m.get("Position", "")
            mutasi = m.get("Mutation", "")
            kelas_obat = m.get("DrugClass", "")
            
            raw_interpretasi = m.get("Interpretasi Lengkap", m.get("Interpretation", "-")) or "-"
            # Split per titik koma, lalu strip spasi
            interpretasi_parts = [p.strip() for p in raw_interpretasi.split(";") if p.strip()]
            
            sample_id = m.get("SampleID", "-")
            reference = m.get("Reference", "-")
            
            # Tulis field standard:
            text_object.textLine(f"Position    : {posisi}")
            text_object.textLine(f"Mutation    : {mutasi}")
            text_object.textLine(f"Drug Class  : {kelas_obat}")
            
            # Tulis interpretasi per‐baris:
            text_object.textLine("Interpretasi:")
            for part in interpretasi_parts:
                text_object.textLine(f"    {part}")
            
            # Lanjutkan dengan Sample ID dan Reference
            text_object.textLine(f"Sample ID   : {sample_id}")
            text_object.textLine(f"Reference   : {reference}")
            
            # Baris kosong pemisah, kecuali ini entry terakhir
            if idx < len(muts) - 1:
                text_object.textLine("")

    c.drawText(text_object)
    c.showPage()
    c.save()
    buf.seek(0)
    return buf