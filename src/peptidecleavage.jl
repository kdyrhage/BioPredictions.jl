"""
    cleave(seq::AminoAcidSequence)

Predict cleavage sites. Returns a vector of cleavage site indices (excluding gaps) and a vector containing the resulting peptides. Currently only supports trypsin.
"""
function cleave(seq::LongAminoAcidSeq)
    seq = ungap(seq)
    seq[end] != AA_Term && push!(seq, AA_Term)
    peptides = LongAminoAcidSeq[]
    cleavagesites = Int[]
    last = 0
    for i in 1:length(seq)
        # window = seq[i-1:i+1]
        if trypsin(seq, i)
            push!(peptides, seq[last+1:i])
            push!(cleavagesites, i)
            last = i
        end
    end
    last < lastindex(seq) && push!(peptides, seq[last:end-1])
    return (positions = cleavagesites, peptides = peptides)
end

function trypsin(seq, i)
    if i-1 <= 0 || i+1 > lastindex(seq)
        return false
    end
    window = seq[i-1:i+1]
    if window == aa"MRP" || window == aa"WKP"
        return true
    elseif window == aa"CKD" || window == aa"DKD" ||
           window == aa"CKH" || window == aa"CKY" ||
           window == aa"CRK" || window == aa"RRR" ||
           window == aa"RRH"
        return false
    elseif (window[2] == AA_K || window[2] == AA_R) && window[3] != AA_P
        return true
    end
    return false
end
