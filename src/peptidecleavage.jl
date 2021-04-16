"""
    cleave(seq::AminoAcidSequence, peptidase::Function = trypsin)

Predict cleavage sites. Returns a vector of cleavage site indices (excluding gaps) and a vector containing the resulting peptides. `protease` is a function that takes two arguments (`seq` and an index `i`) and returns `true` if the sequence would be cleaved between position `i` and `i+1` by the given peptidase.
Available (but unexported) options for `protease` are
- `trypsin`
- `proteinaseK`
- `lysC`
- `proline_endopeptidase`
For details on cleavage sites, see https://web.expasy.org/peptide_cutter/peptidecutter_enzymes.html
"""
function cleave(seq::LongAminoAcidSeq, protease::Function = trypsin)
    seq = ungap(seq)
    cleavagesites = Int[]
    last = 0
    for i in 1:length(seq)
        if protease(seq, i)
            push!(cleavagesites, i)
            last = i
        end
    end
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
    i >= lastindex(seq) && return false
    if seq[i] in [AA_A, AA_E, AA_F, AA_I, AA_L, AA_J, AA_T, AA_V, AA_W, AA_Y]
        return true
    end
    return false
end

function lysC(seq, i)
    i >= lastindex(seq) && return false
    seq[i] == AA_K && return true
    return false
end

function proline_endopeptidase(seq, i)
    i >= lastindex(seq) && return false
    if seq[i] == AA_P && seq[i-1] in [AA_H, AA_K, AA_R] && seq[i+1] != AA_P
        return true
    end
    return false
end
