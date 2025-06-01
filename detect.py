from align import align_sequences
from translate import translate_if_nucleotide
from database import get_db, MutationRule

def get_mutation_rules():
    db = next(get_db())
    rules = db.query(MutationRule).all()
    rule_dict = {(r.position, r.ref_aa, r.alt_aa): {"drug": r.drug_class, "description": r.description} for r in rules}
    return rule_dict

def detect_mutations(ref_prot: str, query_prot: str, mode="global", mutation_rules=None) -> list[dict]:
    aligned_ref, aligned_query = align_sequences(ref_prot, query_prot, mode)
    muts = []
    ref_pos = 0
    for r, q in zip(aligned_ref, aligned_query):
        if r != "-":
            ref_pos += 1
            rule = mutation_rules.get((ref_pos, r, q))
            if rule:
                muts.append({
                    "Position": ref_pos,
                    "Mutation": f"{r}{ref_pos}{q}",
                    "DrugClass": rule["drug"],
                    "Interpretation": rule["description"]
                })
    return muts