import parasail

def align_sequences(ref: str, query: str, mode: str = "global") -> tuple[str,str]:
    if mode == "global":
        aln = parasail.nw_trace_striped_16(query, ref, 10, 1, parasail.blosum62)
    else:
        aln = parasail.sw_trace_striped_16(query, ref, 10, 1, parasail.blosum62)
    tb = aln.traceback
    return tb.reference, tb.query