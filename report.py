import pandas as pd, io
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter

def to_csv(muts: list[dict]) -> bytes:
    return pd.DataFrame(muts).to_csv(index=False).encode()

def to_pdf(muts: list[dict]) -> io.BytesIO:
    """
    Buat PDF dengan format per-entry seperti:

    Position    : ...
    Mutation    : ...
    Drug Class  : ...
    Interpretasi: ...
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
            # Ambil masing‐masing field, dengan fallback ke "-" jika kosong atau None
            posisi = m.get("Position", "")
            mutasi = m.get("Mutation", "")
            kelas_obat = m.get("DrugClass", "")
            
            # Ambil interpretasi lengkap jika ada, jika tidak ambil Interpretation
            interpretasi = m.get("Interpretasi Lengkap", m.get("Interpretation", "-"))
            if not interpretasi:
                interpretasi = "-"
            
            sample_id = m.get("SampleID", "-")
            reference = m.get("Reference", "-")
            
            # Tulis satu baris per field
            text_object.textLine(f"Position    : {posisi}")
            text_object.textLine(f"Mutation    : {mutasi}")
            text_object.textLine(f"Drug Class  : {kelas_obat}")
            text_object.textLine(f"Interpretasi: {interpretasi}")
            text_object.textLine(f"Sample ID   : {sample_id}")
            text_object.textLine(f"Reference   : {reference}")
            
            # Jika masih ada entry berikutnya, beri satu baris kosong sebagai pemisah
            if idx < len(muts) - 1:
                text_object.textLine("")  # baris kosong

    c.drawText(text_object)
    c.showPage()
    c.save()
    buf.seek(0)
    return buf