import parasail
from Bio import pairwise2
from Bio.Align import substitution_matrices 

def align_sequences(ref: str, query: str, mode: str = "global", method: str = "parasail") -> tuple[str, str]:
    if not ref or not query:
        if mode == "global":
            aligned_ref_on_empty = ref + '-' * len(query)
            aligned_query_on_empty = '-' * len(ref) + query
            return aligned_ref_on_empty, aligned_query_on_empty
        else: 
            return "", ""

    if method == "parasail":
        gap_open_parasail = 10
        gap_extend_parasail = 1
        
        matrix_parasail = parasail.blosum62 

        if mode == "global":
            aln = parasail.nw_trace_striped_16(query, ref, gap_open_parasail, gap_extend_parasail, matrix_parasail)
        elif mode == "local":
            aln = parasail.sw_trace_striped_16(query, ref, gap_open_parasail, gap_extend_parasail, matrix_parasail)
        else:
            raise ValueError(f"Mode alignment tidak didukung: {mode} untuk parasail. Pilih 'global' atau 'local'.")
        
        if aln.traceback is None:
            if mode == "global":
                return (ref + '-' * len(query), query + '-' * len(ref))
            else:
                return ("", "")

        return str(aln.traceback.reference), str(aln.traceback.query)

    elif method == "biopython":
        try:
            blosum62 = substitution_matrices.load("BLOSUM62")
        except FileNotFoundError:
            raise RuntimeError("Matriks BLOSUM62 tidak ditemukan oleh Biopython. Pastikan Biopython terinstal dengan benar beserta file datanya.")

        gap_open_biopython = -10
        gap_extend_biopython = -1

        if mode == "global":
            alignments = pairwise2.align.globalds(query, ref, blosum62, gap_open_biopython, gap_extend_biopython)
        elif mode == "local":
            alignments = pairwise2.align.localds(query, ref, blosum62, gap_open_biopython, gap_extend_biopython)
        else:
            raise ValueError(f"Mode alignment tidak didukung: {mode} untuk biopython. Pilih 'global' atau 'local'.")

        if not alignments:
            if mode == "global":
                return (ref + '-' * len(query), query + '-' * len(ref)) 
            else: # local
                return ("", "") 

        aligned_query = alignments[0][0]
        aligned_ref = alignments[0][1]
        
        return aligned_ref, aligned_query
    else:
        raise ValueError(f"Metode alignment tidak didukung: {method}. Pilih 'parasail' atau 'biopython'.")