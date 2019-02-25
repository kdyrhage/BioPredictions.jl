module ProteinStats

export pI, weight, gravy, cleave

using BioSequences
using StatsBase

include("tables.jl")
include("pI.jl")
include("molecularweight.jl")
include("gravy.jl")
include("peptidecleavage.jl")

end
