from align import align_sequences
from translate import translate_if_nucleotide
from database import get_db, MutationRule

def get_mutation_rules():
    db = next(get_db())
    rules = db.query(MutationRule).all()
    rule_dict = {}
    for r in rules:
        rule_dict[(r.position, r.ref_aa.upper(), r.alt_aa.upper())] = {
            "drug": r.drug_class,
            "description": r.description
        }
    return rule_dict

def detect_mutations(ref_prot: str, query_prot: str, mode: str = "global", method: str = "parasail", mutation_rules=None) -> list[dict]:
    if mutation_rules is None:
        mutation_rules = {}

    print(f"Number of mutation rules loaded: {len(mutation_rules)}")
    if mutation_rules:
        print("Sample mutation rules:")
        for i, (key, value) in enumerate(list(mutation_rules.items())[:5]):
            print(f"  {key}: {value}")

    aligned_ref, aligned_query = align_sequences(ref_prot.upper(), query_prot.upper(), mode, method=method)
    
    muts = []
    all_differences = []
    ref_pos = 0
    
    for i, (r_aa, q_aa) in enumerate(zip(aligned_ref, aligned_query)):
        if r_aa != "-":  
            ref_pos += 1

        if r_aa != "-" and q_aa != "-" and r_aa != q_aa:
            all_differences.append({
                "position": ref_pos,
                "reference": r_aa,
                "query": q_aa,
                "mutation": f"{r_aa.upper()}{ref_pos}{q_aa.upper()}"
            })
            
            rule_key = (ref_pos, r_aa.upper(), q_aa.upper())
            rule = mutation_rules.get(rule_key)
            
            if rule:
                muts.append({
                    "Position": ref_pos,
                    "Mutation": f"{r_aa.upper()}{ref_pos}{q_aa.upper()}",
                    "DrugClass": rule["drug"],
                    "Interpretation": rule["description"]
                })
            else:
                print(f"No rule found for mutation: {rule_key}")
    
    print(f"\nTotal differences found: {len(all_differences)}")
    print("All differences:")
    for diff in all_differences[:20]:  # Show first 20 differences
        print(f"  Position {diff['position']}: {diff['reference']} -> {diff['query']} ({diff['mutation']})")
    
    if len(all_differences) > 20:
        print(f"  ... and {len(all_differences) - 20} more differences")
    
    print(f"\nMutations with rules found: {len(muts)}")
    
    return muts