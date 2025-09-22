
function discretize(d::SAX, series::Vector{Float64})
	n = length(series)
	# todo: make a `normalize()` fn
	if d.normalize
		μ = mean(series)
		σ = std(series, corrected=false)
		# handle constant series
		if σ < 1e-10
			middle_symbol = Char('a' + div(d.alphabet_size, 2))
			return fill(middle_symbol, d.nsegments)
		end
		normalized = (series .- μ) ./ σ
	else
		normalized = series
	end
	paa_values = paa(normalized, d.nsegments)
	symbols = [value_to_symbol(paa_values[i], d.breakpoints) for i in 1:d.nsegments]
	return symbols
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

function discretize(d::DeltaDiscretizer, series::Vector{Float64})
	diff_series = series
	for _ in 1:d.order
		diff_series = diff(diff_series)
	end
	isnothing(d.base_discretizer) || return discretize(d.base_discretizer, diff_series)
	isempty(diff_series) && return Char[]
	min_val, max_val = extrema(diff_series)
	range_val = max_val - min_val
	if range_val < 1e-10
		middle_symbol = Char('a' + div(d.alphabet_size, 2))
		return fill(middle_symbol, length(diff_series))
	end
	symbols = Vector{Char}(undef, length(diff_series))
	for i in 1:length(diff_series)
		normalized = (diff_series[i] - min_val) / range_val
		bin_idx = min(floor(Int, normalized * d.alphabet_size), d.alphabet_size - 1)
		symbols[i] = Char('a' + bin_idx)
	end
	return symbols
end

