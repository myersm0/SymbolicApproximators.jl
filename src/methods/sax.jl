
struct SAX <: SymbolicApproximator
	w::Union{Nothing, Int}  # word size
	a::UInt32               # alphabet size
	β::Vector{Float64}      # breakpoints
end

# todo: enforce cardinality >= 2
function SAX(w::Integer, a::Integer)
	β = quantile.(Normal(), (1:a-1) ./ a)
	return SAX(w, a, β)
end

function _encode_segment(sax::SAX, values::AbstractVector{Float64})
	return _encode_segment(PAA(), values)
end

WordStyle(::SAX) = SimpleWord()

