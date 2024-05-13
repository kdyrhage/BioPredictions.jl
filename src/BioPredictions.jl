module BioPredictions

export pI, weight, gravy, cleave, cai

using BioSequences
using StatsBase
using GenomicAnnotations
using Kmers

include("tables.jl")
include("pI.jl")
include("molecularweight.jl")
include("gravy.jl")
include("peptidecleavage.jl")
include("codonadaptationindex.jl")

end
