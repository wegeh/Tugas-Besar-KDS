import streamlit as st
from Bio import SeqIO
import pandas as pd
from translate import translate_if_nucleotide
from detect import detect_mutations, get_mutation_rules
from report import to_csv, to_pdf
from database import init_db, get_db, Sample, Mutation, ReferenceSequence

init_db()

st.title("Deteksi Mutasi Resistensi Obat pada Gen pol HIV-1")
mutation_rules = get_mutation_rules()
uploaded_file = st.file_uploader("Unggah file FASTA (sample pasien)", type=["fasta", "fa"])
mode = st.radio("Mode alignment", options=["global", "local"], index=0)

if st.button("Analyze"):
    if not uploaded_file:
        st.warning("Harap unggah file FASTA sample pasien!")
    else:
        try:
            db = next(get_db())
            all_muts = []
            sample_ids_seen = set()

            # Ambil semua referensi dari database internal
            reference_db_seqs = db.query(ReferenceSequence).all()
            if not reference_db_seqs:
                st.error("Database referensi kosong. Import dulu file FASTA referensi ke DB.")
                st.stop()

            # Loop setiap pasien
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

                # Untuk setiap referensi di database, lakukan komparasi
                for ref_db in reference_db_seqs:
                    ref_prot = ref_db.sequence  # sudah protein, atau translate jika perlu
                    muts = detect_mutations(ref_prot, query_prot, mode, mutation_rules=mutation_rules)
                    for m in muts:
                        m["SampleID"] = rec.id
                        m["Reference"] = ref_db.name
                        db_mut = Mutation(
                            sample_id=samp.id,
                            reference_id=ref_db.id,
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
                st.success(f"Ditemukan {len(all_muts)} mutasi pada {len(sample_ids_seen)} sampel.")
                st.dataframe(df)
                st.download_button("Download CSV", to_csv(all_muts), "mutations.csv")
                st.download_button("Download PDF", to_pdf(all_muts), "mutations.pdf")
            else:
                st.info("Tidak ada mutasi resistensi terdeteksi pada seluruh sampel.")

        except Exception as e:
            st.error(f"Terjadi error saat analisis: {str(e)}")
