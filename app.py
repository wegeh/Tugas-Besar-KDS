import streamlit as st
from Bio import SeqIO
import pandas as pd
import io
from translate import translate_if_nucleotide
from detect import detect_mutations, get_mutation_rules
from report import to_csv, to_pdf
from database import init_db, get_db, Sample, Mutation, ReferenceSequence

DRUG_FULL_NAMES = {
    "ABC": "Abacavir", "AZT": "Zidovudine", "FTC": "Emtricitabine", "3TC": "Lamivudine", "TDF": "Tenofovir DF",
    "DDI": "Didanosine", "D4T": "Stavudine",
    "DOR": "Doravirine", "EFV": "Efavirenz", "ETR": "Etravirine", "NVP": "Nevirapine", "RPV": "Rilpivirine",
    "DLV": "Delavirdine",
    "ATV": "Atazanavir", "ATV/R": "Atazanavir/ritonavir",
    "DRV": "Darunavir", "DRV/R": "Darunavir/ritonavir",
    "LPV": "Lopinavir", "LPV/R": "Lopinavir/ritonavir",
    "FPV": "Fosamprenavir", "IDV": "Indinavir", "NFV": "Nelfinavir", "SQV": "Saquinavir", "TPV": "Tipranavir",
    "BIC": "Bictegravir", "CAB": "Cabotegravir", "DTG": "Dolutegravir", "EVG": "Elvitegravir", "RAL": "Raltegravir",
    "LEN": "Lenacapavir",
}

def get_qualitative_resistance_level(score_str: str) -> str:
    """Mengubah skor numerik (string) menjadi deskripsi kualitatif."""
    try:
        score = int(score_str)
        if score < 0:
            return "    "
        elif 0 <= score <= 9:
            return "Suseptibel"
        elif 10 <= score <= 14:
            return "Potensi Resistensi Rendah"
        elif 15 <= score <= 29:
            return "Resistensi Rendah"
        elif 30 <= score <= 59:
            return "Resistensi Menengah"
        elif score >= 60:
            return "Resistensi Tinggi"
        else:
            return f"Skor {score}"
    except ValueError:
        return score_str 
    except TypeError:
        return "N/A"

def format_interpretation_for_display(interpretation_scores_str: str) -> str:
    """Mengubah string skor menjadi deskripsi kualitatif yang lebih ramah pengguna."""
    if not interpretation_scores_str or interpretation_scores_str == "-":
        return "-"

    parts = interpretation_scores_str.split(',')
    qualitative_descriptions = []
    for part in parts:
        try:
            drug_code_raw, score_value_raw = part.split(':', 1)
            drug_code = drug_code_raw.strip().upper() # Untuk lookup di DRUG_FULL_NAMES
            score_value = score_value_raw.strip()
            
            # Cari nama lengkap, fallback ke kode jika tidak ada
            full_drug_name = DRUG_FULL_NAMES.get(drug_code, drug_code_raw.strip()) 
            
            qualitative_level = get_qualitative_resistance_level(score_value)
            qualitative_descriptions.append(f"{full_drug_name}: {qualitative_level} (Skor: {score_value})")
        except ValueError:
            qualitative_descriptions.append(part.strip()) # Tampilkan apa adanya jika format salah

    return "; ".join(qualitative_descriptions)

init_db()
db_session = next(get_db())

st.markdown("<h1 style='text-align: center; margin-bottom: 40px;'>Deteksi Mutasi Resistensi Obat pada Gen pol HIV-1</h1>", unsafe_allow_html=True)

mutation_rules = get_mutation_rules()

try:
    all_references_from_db = db_session.query(ReferenceSequence).all()
    reference_names = [ref.name for ref in all_references_from_db]
except Exception as e:
    st.error(f"Gagal memuat daftar referensi dari database: {e}")
    reference_names = []

uploaded_file = st.file_uploader("Unggah file FASTA (sample pasien)", type=["fasta", "fa"])
mode = st.radio("Mode alignment (Global/Lokal):", options=["global", "local"], index=0, horizontal=True)

alignment_method_choice = st.radio(
    "Pilih Pustaka Alignment:",
    options=["parasail", "biopython"],
    index=0,
    horizontal=True,
    help="Pilih 'parasail' untuk kecepatan (default) atau 'biopython' untuk alternatif."
)

selected_reference_name = None
if reference_names:
    selected_reference_name = st.selectbox("Pilih Sekuens Referensi:", options=reference_names)
else:
    st.warning("Tidak ada sekuens referensi di database. Mohon impor terlebih dahulu.")


if st.button("Analyze"):
    if not uploaded_file:
        st.warning("Harap unggah file FASTA sample pasien!")
    elif not selected_reference_name and reference_names:
        st.warning("Harap pilih sekuens referensi!")
    elif not reference_names:
        st.error("Database referensi kosong. Tidak dapat melanjutkan analisis.")
    else:
        try:
            db = db_session 

            all_muts = []
            sample_ids_seen = set()

            chosen_reference_db = db.query(ReferenceSequence).filter_by(name=selected_reference_name).first()
            if not chosen_reference_db:
                st.error(f"Referensi '{selected_reference_name}' tidak ditemukan di database.")
                st.stop()

            uploaded_file.seek(0)
            with io.TextIOWrapper(uploaded_file, encoding="utf-8-sig") as fasta_textio:
                for rec in SeqIO.parse(fasta_textio, "fasta"):
                    if rec.id in sample_ids_seen:
                        continue
                    sample_ids_seen.add(rec.id)
                    samp = db.query(Sample).filter_by(sample_id=rec.id).first()
                    if not samp:
                        samp = Sample(
                            sample_id=rec.id,
                            description=rec.description or "",
                            sequence=str(rec.seq)
                        )
                        db.add(samp)
                        db.commit()
                        db.refresh(samp)

                    query_prot = translate_if_nucleotide(rec)
                    ref_prot = chosen_reference_db.sequence

                    muts = detect_mutations(
                        ref_prot,
                        query_prot,
                        mode,
                        method=alignment_method_choice,
                        mutation_rules=mutation_rules
                    )
                    for m in muts:
                        m["SampleID"] = rec.id
                        m["Reference"] = chosen_reference_db.name
                        db_mut = Mutation(
                            sample_id=samp.id,
                            reference_id=chosen_reference_db.id,
                            position=m.get("Position", 0),
                            mutation=m.get("Mutation", ""),
                            drug_class=m.get("DrugClass", "unknown"),
                            interpretation=m.get("Interpretation", "-")
                        )
                        db.add(db_mut)
                    all_muts.extend(muts)

                db.commit()

                if all_muts:
                    # 1) Buat DataFrame dari list semua mutasi
                    df = pd.DataFrame(all_muts)

                    # 2) Tambahkan kolom baru dengan versi "dijelaskan" dari Interpretation
                    df["Interpretasi Lengkap"] = df["Interpretation"].apply(
                        lambda raw: format_interpretation_for_display(raw)
                    )

                    # 3) Tampilkan kolom-kolom utama, termasuk kolom yang sudah dijelaskan
                    #    Anda bisa memilih kolom apa saja yang ingin ditampilkan di tabel Streamlit
                    col_order = [
                        "SampleID", "Reference", "Position", "Mutation", "DrugClass",
                        "Interpretasi Lengkap"
                    ]
                    st.success(
                        f"Ditemukan {len(all_muts)} mutasi pada {len(sample_ids_seen)} sampel "
                        f"terhadap referensi '{selected_reference_name}' menggunakan metode '{alignment_method_choice}'."
                    )
                    st.dataframe(df[col_order])

                    # 4) Tombol download CSV / PDF tetap menggunakan `all_muts`, 
                    #    tetapi kalau mau, Anda bisa menyertakan kolom "Interpretasi Lengkap" di file-nya.
                    col1, col2, col3, col4, col5 = st.columns(5)
                    with col1:
                        st.download_button(
                            "Download CSV",
                            to_csv(df.to_dict("records")),  # convert DataFrame ke list of dict
                            f"mutations_{selected_reference_name}_{alignment_method_choice}.csv",
                            key="csv_download"
                        )
                    with col2:
                        st.download_button(
                            "Download PDF",
                            to_pdf(df.to_dict("records")),
                            f"mutations_{selected_reference_name}_{alignment_method_choice}.pdf",
                            key="pdf_download"
                        )
                else:
                    st.info(
                        f"Tidak ada mutasi resistensi terdeteksi pada seluruh sampel terhadap referensi "
                        f"'{selected_reference_name}' menggunakan metode '{alignment_method_choice}'."
                    )

        except Exception as e:
            st.error(f"Terjadi error saat analisis: {str(e)}")
        finally:
            if db_session:
                db_session.close()