"""
    cleave(seq::AminoAcidSequence, peptidase::Function = trypsin)

Predict cleavage sites. Returns a vector of cleavage site indices (excluding gaps) and a vector containing the resulting peptides. `protease` is a function that takes two arguments (`seq` and an index `i`) and returns `true` if the sequence would be cleaved at position `i` by the given peptidase.
Available options for `protease` are
    - `trypsin`
    - `proteinaseK`
"""
function cleave(seq::LongAminoAcidSeq, protease::Function = trypsin)
    seq = ungap(seq)
    # seq[end] != AA_Term && push!(seq, AA_Term)
    # peptides = LongAminoAcidSeq[]
    cleavagesites = Int[]
    last = 0
    for i in 1:length(seq)
        if protease(seq, i)
            # push!(peptides, seq[last+1:i])
            push!(cleavagesites, i)
            last = i
        end
    end
    # last < lastindex(seq) && push!(peptides, seq[last+1:end])
    peptides = getpeptides(seq, cleavagesites)
    return (positions = cleavagesites, peptides = peptides)
end


function getpeptides(seq, sites)
    peptides = Array{LongAminoAcidSeq}(undef, length(sites)+1)
    last = 1
    for (i, site) in enumerate(sites)
        peptides[i] = seq[last:site]
        last = site + 1
    end
    peptides[end] = seq[last:end]
    return peptides
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

function proteinaseK(seq, i)
    if i >= lastindex(seq) return false
    end
    if seq[i] in [AA_A, AA_E, AA_F, AA_I, AA_L, AA_J, AA_T, AA_V, AA_W, AA_Y]
        return true
    end
    return false
end
