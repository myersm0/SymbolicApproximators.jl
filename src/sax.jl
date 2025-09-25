
struct SAX{T, A <: AbstractVector{T}} <: SymbolicApproximator{T, A}
	w::Int              # word size
	α::A                # alphabet
	β::Vector{Float64}  # breakpoints
end

function SAX(w::Integer, α::AbstractVector{T}) where T
	cardinality = length(α)
	cardinality > 1 || error("cardinality must be at least 2")
	β = quantile.(Normal(), (1:cardinality-1) ./ cardinality)
	return SAX{T, typeof(α)}(w, α, β)
end

function SAX(w::Integer, cardinality::Integer)
	# todo: more helpful error msg
	1 < cardinality <= 26 || error("invalid alphabet size")
	α = 'a':('a' + cardinality - 1)
	return SAX(w, α)
end

function _encode_segment(disc::SAX, values::Vector{Float64})
	μ = mean(values)
	σ = std(values, corrected=false)
	# handle constant values
	if σ < 1e-10
		middle_symbol = Char('a' + div(length(disc.α), 2))
		return fill(middle_symbol, d.w)
	end
	normalized = (values .- μ) ./ σ
	paa_values = paa(normalized, disc.w)
	return Word(disc, paa_values)
end

