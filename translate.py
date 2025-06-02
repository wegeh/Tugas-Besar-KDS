
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def translate_if_nucleotide(rec: SeqRecord) -> str:
    seq = str(rec.seq).upper()
    if set(seq).issubset(set("ACGTU")) and len(seq) > 20:  
        best_translation = ""
        best_length = 0

        for frame in range(3):
            try:
                translated = str(Seq(seq[frame:]).translate())
                translated = translated.rstrip('*')
                if len(translated) > best_length:
                    best_translation = translated
                    best_length = len(translated)
                    
            except Exception as e:
                print(f"Translation error in frame {frame}: {e}")
                continue
        
        if best_translation:
            print(f"Translated sequence length: {len(best_translation)}")
            print(f"First 50 amino acids: {best_translation[:50]}")
            return best_translation
        else:
            print("Translation failed, returning original sequence")
            return seq
    
    print(f"Sequence appears to be protein, length: {len(seq)}")
    return seq