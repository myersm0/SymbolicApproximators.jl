
abstract type SymbolicDiscretizer{T, A} end

eltype(::SymbolicDiscretizer{T, A}) where {T, A} = T

alphabet(disc::SymbolicDiscretizer) where {T, A} = disc.α
breakpoints(disc::SymbolicDiscretizer) where {T, A} = disc.β
cardinality(disc::SymbolicDiscretizer) where {T, A} = length(alphabet(disc))

function discretize(::SymbolicDiscretizer, values::AbstactVector) end




