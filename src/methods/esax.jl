
struct ESAX <: SymbolicApproximator
	w::Union{Nothing, Int}  # word size
	a::UInt32               # alphabet size
	β::Vector{Float64}      # breakpoints
end

function ESAX(w::Integer, a::Integer)
	β = [-Inf; quantile.(Normal(), (1:a-1) ./ a)]
	return ESAX(w, a, β)
end

function _encode_segment(sa::ESAX, values::AbstractVector{<:Real})
	pmin = argmin(values)
	pmax = argmax(values)
	meanval = mean(values)
	pmid = (1 + length(values)) ÷ 2   # midpoint of segment (ESAX Eq. 2)
	perm = sortperm([pmin, pmid, pmax])
	return SVector{3, Float64}([values[pmin], meanval, values[pmax]][perm])
end

WordStyle(::ESAX) = MultiWord{3}()



