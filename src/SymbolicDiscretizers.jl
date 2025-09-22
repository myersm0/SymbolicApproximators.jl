
module SymbolicDiscretizers

using LinearAlgebra
using Distributions
using StatsBase

abstract type SymbolicDiscretizer end
export SymbolicDiscretizer

export discretize, reconstruct, distance

include("utilities.jl")
export build_distance_table, numerosity_reduction

# main algorithms
include("sax.jl")
include("ordinal.jl")
include("delta.jl")
export SAX, OrdinalDiscretizer, DeltaDiscretizer

include("windows.jl")
export StreamingDiscretizer, update!, get_symbols

end
