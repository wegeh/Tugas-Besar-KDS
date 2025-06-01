from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def translate_if_nucleotide(rec: SeqRecord) -> str:
    seq = str(rec.seq).upper()
    if set(seq).issubset(set("ACGTU")):
        return str(Seq(seq).translate(to_stop=True))
    return seq
