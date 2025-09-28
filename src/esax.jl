
struct ESAX{T} <: SymbolicApproximator{T}
	w::Union{Nothing, Int}  # word size
	α::Vector{T}            # alphabet
	β::Vector{Float64}      # breakpoints
end

function ESAX(w::Integer, α::AbstractVector{T}) where T
	cardinality = length(α)
	β = [-Inf; quantile.(Normal(), (1:cardinality-1) ./ cardinality)]
	return ESAX{T, typeof(α)}(w, α, β)
end

function ESAX(w::Integer, cardinality::Integer)
	α = 'a':('a' + cardinality - 1)
	return ESAX(w, α)
end

function _encode_segment(sa::ESAX, values::AbstractVector{Float64})
	pmin = argmin(values)
	pmax = argmax(values)
	meanval = mean(values)
	pmid = (1 + length(values)) ÷ 2   # midpoint of segment (ESAX Eq. 2)
	perm = sortperm([pmin, pmid, pmax])
	return SVector{3, Float64}([values[pmin], meanval, values[pmax]][perm])
end

function Word(sa::SA, values::AbstractVector) where SA <: ESAX
	β = breakpoints(sa)
	n = length(values)
	indices = Vector{SVector{3, Int}}(undef, n)
	for (i, v) in enumerate(values)
		indices[i] = SVector{3, Int}(searchsortedlast(β, vi) for vi in v)
	end
	return Word{SA, SVector{3, Int}, 3}(Ref(sa), indices, n)
end

width(::ESAX) = 3






