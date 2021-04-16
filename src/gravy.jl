"""
    gravy(sequence)

Calculate GRAVY index of `sequence`.
Source: https://www.ncbi.nlm.nih.gov/pubmed/7108955
"""
function gravy(sequence)
    aacomposition = BioSequences.composition(sequence)
    gravy = 0.0
    for (aa, n) in aacomposition
        gravy += n * get(hydropathy, aa, 0.0)
    end
    return gravy / length(sequence)
end
