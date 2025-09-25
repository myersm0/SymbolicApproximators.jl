
struct PAA{T <: Real} <: ContinuousApproximator{T}
	w::Union{Nothing, Int}  # word size
end

PAA() = PAA{Float64}(nothing)
PAA(w::Integer) = PAA{Float64}(w)

function _encode_segment(::PAA, values::AbstractVector{<:Real})
	return mean(values)
end

