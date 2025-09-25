
abstract type SymbolicApproximator{T, A} end

Base.eltype(::SymbolicApproximator{T, A}) where {T, A} = T
alphabet(disc::SymbolicApproximator) where {T, A} = disc.α
breakpoints(disc::SymbolicApproximator) where {T, A} = disc.β
cardinality(disc::SymbolicApproximator) where {T, A} = length(alphabet(disc))

function discretize(::SymbolicApproximator, values::AbstactVector) end




