module ProteinStats

export pI

using BioSequences
using StatsBase

include("tables.jl")
include("pI.jl")
include("molecularweight.jl")
include("gravy.jl")

end
