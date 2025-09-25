
struct SAX{T, A <: AbstractVector{T}} <: SymbolicDiscretizer{T, A}
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

function discretize(disc::SAX, series::Vector{Float64})
	μ = mean(series)
	σ = std(series, corrected=false)
	# handle constant series
	if σ < 1e-10
		middle_symbol = Char('a' + div(length(disc.α), 2))
		return fill(middle_symbol, d.w)
	end
	normalized = (series .- μ) ./ σ
	paa_values = paa(normalized, disc.w)
#	symbols = [value_to_symbol(paa_values[i], disc.breakpoints) for i in 1:disc.w]
	word = Word(disc, paa_values)
	return word
end

function distance(disc::SAX, symbols1::Vector{Char}, symbols2::Vector{Char})
	length(symbols1) == length(symbols2) || error("Symbol sequences must have same length")
	dist_sum = 0.0
	for (a, b) in zip(symbols1, symbols2)
		dist_sum += distance(a, b, disc.α)^2
	end
	return sqrt(dist_sum)
end

function reconstruct(sax::SAX, symbols::Vector{Char}, original_length::Int)
	paa_values = zeros(length(symbols))
	for i in 1:length(symbols)
		idx = Int(symbols[i] - 'a')
		if idx == 0
			paa_values[i] = sax.breakpoints[1] - 0.5
		elseif idx >= length(sax.breakpoints)
			paa_values[i] = sax.breakpoints[end] + 0.5
		else
			paa_values[i] = (sax.breakpoints[idx] + sax.breakpoints[idx+1]) / 2
		end
	end
	segment_length = original_length / length(symbols)
	reconstructed = Float64[]
	for i in 1:length(symbols)
		repeat_count = round(Int, segment_length)
		if i == length(symbols)
			repeat_count = original_length - length(reconstructed)
		end
		append!(reconstructed, fill(paa_values[i], repeat_count))
	end
	return reconstructed[1:original_length]
end

