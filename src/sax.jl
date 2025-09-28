
struct SAX{T} <: SymbolicApproximator{T}
	w::Union{Nothing, Int}  # word size
	α::Vector{T}            # alphabet
	β::Vector{Float64}      # breakpoints
end

# todo: enforce cardinality >= 2
function SAX(w::Integer, α::AbstractVector{T}) where T
	cardinality = length(α)
	β = [-Inf; quantile.(Normal(), (1:cardinality-1) ./ cardinality)]
	return SAX{T}(w, α, β)
end

# todo: enforce alphabet size <= 26 for this method
function SAX(w::Integer, cardinality::Integer)
	α = 'a':('a' + cardinality - 1)
	return SAX(w, α)
end

function _encode_segment(sax::SAX, values::AbstractVector{Float64})
	return _encode_segment(PAA(), values)
end

WordStyle(::SAX) = SimpleWord()

