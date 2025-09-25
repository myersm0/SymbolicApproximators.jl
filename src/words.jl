
struct Word{D, T}
	disc::Ref{D}
	symbols::Vector{T}
end

function Word(disc::D, values) where D
	α = disc.α
	β = disc.β
	T = eltype(α)
	n = length(values)
	symbols = Vector{T}(undef, n)
	@inbounds @simd for i in eachindex(values)
		v = Float64(values[i])
		index = searchsortedlast(β, v)
		symbols[i] = α[index + 1]
	end
	Word{D,T}(Ref(disc), symbols)
end







