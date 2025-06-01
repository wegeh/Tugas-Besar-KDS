from database import init_db, get_db, MutationRule
init_db()
db = next(get_db())

rules = [
    # Pos, Ref, Alt, Drug, Description
    (103, "K", "N", "NNRTI", "High-level resistance to efavirenz/nevirapine"),
    (184, "M", "V", "NRTI",  "High-level resistance to lamivudine/emtricitabine"),
    (65,  "K", "R", "NRTI",  "Resistance to tenofovir"),
    # tambahkan 
]

for pos, ref, alt, drug, desc in rules:
    db.add(MutationRule(position=pos, ref_aa=ref, alt_aa=alt, drug_class=drug, description=desc))
db.commit()
