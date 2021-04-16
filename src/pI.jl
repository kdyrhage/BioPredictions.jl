function charge(sequence, pKvalues = pKtable["solomon"], pH = 7.0)
    aa = BioSequences.composition(sequence)

    charge = 0.0

    charge += (1 / (1 + 10^(1 * (pH - pKvalues[1]))))
    charge += (-1 / (1 + 10^(-1 * (pH - pKvalues[2]))))

    charge += aa[AA_C] * (-1 / (1 + 10^(-1 * (pH - pKvalues[3]))))
    charge += aa[AA_D] * (-1 / (1 + 10^(-1 * (pH - pKvalues[4]))))
    charge += aa[AA_E] * (-1 / (1 + 10^(-1 * (pH - pKvalues[5]))))
    charge += aa[AA_Y] * (-1 / (1 + 10^(-1 * (pH - pKvalues[9]))))
    charge += aa[AA_H] * (1 / (1 + 10^(1 * (pH - pKvalues[6]))))
    charge += aa[AA_K] * (1 / (1 + 10^(1 * (pH - pKvalues[7]))))
    charge += aa[AA_R] * (1 / (1 + 10^(1 * (pH - pKvalues[8]))))

    return charge
end


"""
    pI(sequence, [pKvalues; gamma])

Estimate isoelectric point of `sequence`. `pKvalues` can be a vector containing pK values or a string specifying one of the pK-value sources in `pKtable`. See the docstring for `pKtable` for details. Default is "IPC protein".
"""
pI(sequence::LongAminoAcidSeq, pKvalues::AbstractString; gamma = 0.001) = pI(sequence, pKtable[pKvalues]; gamma = gamma)
function pI(sequence::LongAminoAcidSeq, pKvalues = pKtable["IPC protein"]; gamma = 0.001)
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
