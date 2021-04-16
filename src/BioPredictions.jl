module BioPredictions

export pI, weight, gravy, cleave

using BioSequences

include("tables.jl")
include("pI.jl")
include("molecularweight.jl")
include("gravy.jl")
include("peptidecleavage.jl")

end
