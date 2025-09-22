
module SymbolicDiscretizers

using LinearAlgebra
using Distributions
using StatsBase

export SymbolicDiscretizer, SAX, OrdinalDiscretizer, DeltaDiscretizer
export discretize, reconstruct, distance
export StreamingDiscretizer, update!, get_symbols
export build_distance_table, numerosity_reduction

include("utilities.jl")
include("discretizers.jl")
include("windows.jl")
include("discretize.jl")
include("distances.jl")
include("misc.jl")

end
