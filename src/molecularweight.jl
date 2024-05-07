"""
    weight(sequence; ignoreunknown = false)

Estimate molecular weight of `sequence` in Dalton.

If `ignoreunknown` is set to `false` (default), non-specific residues (X: any AA, Z: glutamine or glutamic acid, B: asparagine or aspartic acid, J: leucine or isoleucine) are calculated using the average of all possible residues in that category. Otherwise they are ignored. Only standard residues are used to calculate the average for X.
"""
function weight(sequence; ignoreunknown = false)
    # aacomposition = BioSequences.composition(sequence)
    aacomposition = countmap(sequence)
    mW = 18.01528 # ≈ weight of one water molecule in Dalton
    for (aa, n) in aacomposition
        if aa ∉ [AA_X, AA_Z, AA_B, AA_J] || !ignoreunknown
            mW += get(aminoacidweights, aa, 0.0) * n
        end
    end
    return mW
end
