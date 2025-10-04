
struct PAA <: ContinuousApproximator
	w::Union{Nothing, Int}  # word size
end

PAA() = PAA(nothing)

function _encode_segment(::PAA, values::AbstractVector)
	return mean(values)
end

WordStyle(::PAA) = SimpleWord()

