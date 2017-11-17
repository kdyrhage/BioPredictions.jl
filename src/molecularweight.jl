function weight(sequence)
    aacomposition = BioSequences.composition(sequence)
    mW = 18.01528 # ≈ weight of one water molecule in Dalton
    for (aa, n) in aacomposition
        mW += get(aminoacidweights, aa, 0) * n
    end
    return mW
end
