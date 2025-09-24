
module SymbolicDiscretizers

using LinearAlgebra
using Distributions
using StatsBase

abstract type SymbolicDiscretizer end
export SymbolicDiscretizer

include("utilities.jl")
export numerosity_reduction

# main algorithms
include("sax.jl")
include("ordinal.jl")
include("delta.jl")
export SAX, OrdinalDiscretizer, DeltaDiscretizer
export discretize, reconstruct, distance

include("windows.jl")
export StreamingDiscretizer, update!, get_symbols

end
