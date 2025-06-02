import streamlit as st
from Bio import SeqIO
import pandas as pd
from translate import translate_if_nucleotide
from detect import detect_mutations, get_mutation_rules
from report import to_csv, to_pdf
from database import init_db, get_db, Sample, Mutation, ReferenceSequence

init_db()
db_session = next(get_db())

st.title("Deteksi Mutasi Resistensi Obat pada Gen pol HIV-1")
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

            for rec in SeqIO.parse(uploaded_file, "fasta"):
                if rec.id in sample_ids_seen:
                    continue
                sample_ids_seen.add(rec.id)

                samp = db.query(Sample).filter_by(sample_id=rec.id).first()
                if not samp:
                    samp = Sample(sample_id=rec.id, description=rec.description or "", sequence=str(rec.seq) )
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
                df = pd.DataFrame(all_muts)
                st.success(f"Ditemukan {len(all_muts)} mutasi pada {len(sample_ids_seen)} sampel terhadap referensi '{selected_reference_name}' menggunakan metode '{alignment_method_choice}'.")
                st.dataframe(df)
                st.download_button("Download CSV", to_csv(all_muts), f"mutations_{selected_reference_name}_{alignment_method_choice}.csv", key="csv_download")
                st.download_button("Download PDF", to_pdf(all_muts), f"mutations_{selected_reference_name}_{alignment_method_choice}.pdf", key="pdf_download")
            else:
                st.info(f"Tidak ada mutasi resistensi terdeteksi pada seluruh sampel terhadap referensi '{selected_reference_name}' menggunakan metode '{alignment_method_choice}'.")

        except Exception as e:
            st.error(f"Terjadi error saat analisis: {str(e)}")
        finally:
            if db_session:
                db_session.close()