
struct OrdinalDiscretizer <: SymbolicDiscretizer
	order::Int		    # d in Keller paper, n in Bandt-Pompe
	delay::Int		    # τ (tau) - time delay
	overlapping::Bool  # whether windows overlap
end

function OrdinalDiscretizer(order::Int; delay::Int = 1, overlapping = true)
	order >= 1 || error("order must be at least 1")
	order <= 8 || @warn "order > 8 may be expensive ($((order+1)!) patterns)"
	delay >= 1 || error("delay must be at least 1")
	return OrdinalDiscretizer(order, delay, overlapping)
end

# efficient pattern encoding using inversions (from Keller & Sinn)
# (this is a helper for discretize())
function ordinal_pattern_to_number(series::AbstractVector{Float64}, t::Int, d::Int, τ::Int)
	# compute inversions i_τ^l(t) for l = 1..d
	inversions = zeros(Int, d)
	for l in 1:d
		count = 0
		for r in 0:(l-1)
			if series[t - r*τ] <= series[t - l*τ]
				count += 1
			end
		end
		inversions[l] = count
	end
	# convert to pattern number using factorial number system
	n = 0
	factorial = 1
	for l in 1:d
		factorial *= (l + 1)
		n += inversions[l] * div(factorial, (l + 1))
	end
	return n
end

function discretize(opd::OrdinalDiscretizer, series::Vector{Float64})
	n = length(series)
	d = opd.order
	τ = opd.delay
	min_length = d * τ + 1
	n >= min_length || error("Series too short: need at least $min_length points")
	
	# determine indices for patterns
	if opd.overlapping
		indices = (d*τ + 1):n  # can compute pattern at each valid position
	else
		indices = (d*τ + 1):(d*τ + 1):(n)  # non-overlapping windows
	end
	
	symbols = Vector{Char}(undef, length(indices))
	for (i, t) in enumerate(indices)
		# get pattern number (0 to (d+1)!-1)
		pattern_num = ordinal_pattern_to_number(series, t, d, τ)
		# map to alphabet (todo: may need modulo for large orders)
		symbols[i] = Char('a' + pattern_num)
	end
	return symbols
end

function permutation_entropy(opd::OrdinalDiscretizer, series::Vector{Float64})
	symbols = discretize(opd, series)
	pattern_counts = Dict{Char, Int}()
	for s in symbols
		pattern_counts[s] = get(pattern_counts, s, 0) + 1
	end
	n = length(symbols)
	entropy = 0.0
	for count in values(pattern_counts)
		p = count / n
		p > 0 && (entropy -= p * log2(p))
	end
	return entropy
end

function permutation_entropy_normalized(opd::OrdinalDiscretizer, series::Vector{Float64})
	H = permutation_entropy(opd, series)
	H_max = log2(factorial(opd.order + 1))
	return H / H_max
end

