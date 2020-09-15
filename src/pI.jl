function charge(sequence, pKvalues = pKtable["solomon"], pH = 7.0)
    # aa = countmap(split(sequence, ""))
    aa = BioSequences.composition(sequence)

    charge = 0.0
    charge += (1 / (1 + 10^(1 * (pH - get(pKvalues, "NTerm", 0.0)))))
    charge += (-1 / (1 + 10^(-1 * (pH - get(pKvalues, "CTerm", 0.0)))))

    charge += aa[AA_R] * (1 / (1 + 10^(1 * (pH - get(pKvalues, AA_R, 0.0)))))
    charge += aa[AA_H] * (1 / (1 + 10^(1 * (pH - get(pKvalues, AA_H, 0.0)))))
    charge += aa[AA_K] * (1 / (1 + 10^(1 * (pH - get(pKvalues, AA_K, 0.0)))))
    charge += aa[AA_D] * (-1 / (1 + 10^(-1 * (pH - get(pKvalues, AA_D, 0.0)))))
    charge += aa[AA_E] * (-1 / (1 + 10^(-1 * (pH - get(pKvalues, AA_E, 0.0)))))
    charge += aa[AA_C] * (-1 / (1 + 10^(-1 * (pH - get(pKvalues, AA_C, 0.0)))))
    charge += aa[AA_Y] * (-1 / (1 + 10^(-1 * (pH - get(pKvalues, AA_Y, 0.0)))))

    return charge
end


"""
    pI(sequence, [pKvalues, gamma])
"""
function pI(sequence::LongAminoAcidSeq, pKvalues = pKtable["solomon"]; gamma = 0.001)
    pH = 7.0

    if charge(sequence, pKvalues, pH) < 0
        pHs = 0:gamma:7
        charges = charge.(Ref(sequence), Ref(pKvalues), pHs)
        pH = pHs[findmin(abs.(charges))[2]]
        return pH
    end

    if charge(sequence, pKvalues, pH) > 0
        pHs = 7:gamma:14
        charges = charge.(Ref(sequence), Ref(pKvalues), pHs)
        pH = pHs[findmin(abs.(charges))[2]]
        return pH
    end

    return pH
end
