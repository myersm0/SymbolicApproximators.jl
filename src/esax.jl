
struct ESAX{T, A <: AbstractVector{T}} <: SymbolicApproximator{T, A}
	w::Union{Nothing, Int}  # word size
	α::A                    # alphabet
	β::Vector{Float64}      # breakpoints
end

function ESAX(w::Integer, α::AbstractVector{T}) where T
	cardinality = length(α)
	cardinality > 1 || error("cardinality must be at least 2")
	β = [-Inf; quantile.(Normal(), (1:cardinality-1) ./ cardinality)]
	return ESAX{T, typeof(α)}(w, α, β)
end

function ESAX(w::Integer, cardinality::Integer)
	# todo: more helpful error msg
	1 < cardinality <= 26 || error("invalid alphabet size")
	α = 'a':('a' + cardinality - 1)
	return ESAX(w, α)
end

# 151 ns for vec, 230 for tuple
function _encode_segment(sa::ESAX, values::AbstractVector{Float64})
	pmin = argmin(values)
	pmax = argmax(values)
	meanval = mean(values)
	pmid = (1 + length(values)) ÷ 2   # midpoint of segment (ESAX Eq. 2)
	perm = sortperm([pmin, pmid, pmax])
	return [values[pmin], meanval, values[pmax]][perm]
end






