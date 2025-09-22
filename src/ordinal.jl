
struct OrdinalDiscretizer <: SymbolicDiscretizer
	alphabet_size::Int
	method::Symbol
	function OrdinalDiscretizer(alphabet_size::Int; method::Symbol=:quantile)
		method in [:rank, :quantile] || error("method must be :rank or :quantile")
		new(alphabet_size, method)
	end
end

function discretize(d::OrdinalDiscretizer, series::Vector{Float64})
	n = length(series)
	if d.method == :rank
		ranks = sortperm(sortperm(series))
		symbols = Vector{Char}(undef, n)
		bins = linspace_bins(1, n, d.alphabet_size)
		for i in 1:n
			bin_idx = searchsortedlast(bins, ranks[i])
			bin_idx = max(1, min(bin_idx, d.alphabet_size))
			symbols[i] = Char('a' + bin_idx - 1)
		end
		return symbols
	else  # :quantile
		quantiles = [i / d.alphabet_size for i in 1:(d.alphabet_size-1)]
		breakpoints = [quantile(series, q) for q in quantiles]
		symbols = Vector{Char}(undef, n)
		for i in 1:n
			symbols[i] = value_to_symbol_ordinal(series[i], breakpoints)
		end
		return symbols
	end
end

function distance(d::OrdinalDiscretizer, symbols1::Vector{Char}, symbols2::Vector{Char})
	# simple Hamming distance for ordinal
	length(symbols1) == length(symbols2) || error("Symbol sequences must have same length")
	return sum(symbols1 .!= symbols2)
end

