
function discretize(sax::SAX, series::Vector{Float64})
	n = length(series)
	if sax.normalize
		μ = mean(series)
		σ = std(series, corrected=false)
		# handle constant series
		if σ < 1e-10
			middle_symbol = Char('a' + div(sax.alphabet_size, 2))
			return fill(middle_symbol, sax.nsegments)
		end
		normalized = (series .- μ) ./ σ
	else
		normalized = series
	end
	paa_values = paa(normalized, sax.nsegments)
	symbols = Vector{Char}(undef, sax.nsegments)
	for i in 1:sax.nsegments
		symbols[i] = value_to_symbol(paa_values[i], sax.breakpoints)
	end
	return symbols
end

function discretize(ord::OrdinalDiscretizer, series::Vector{Float64})
	n = length(series)
	if ord.method == :rank
		ranks = sortperm(sortperm(series))
		symbols = Vector{Char}(undef, n)
		bins = linspace_bins(1, n, ord.alphabet_size)
		for i in 1:n
			bin_idx = searchsortedlast(bins, ranks[i])
			bin_idx = max(1, min(bin_idx, ord.alphabet_size))
			symbols[i] = Char('a' + bin_idx - 1)
		end
		return symbols
	else  # :quantile
		quantiles = [i / ord.alphabet_size for i in 1:(ord.alphabet_size-1)]
		breakpoints = [quantile(series, q) for q in quantiles]
		symbols = Vector{Char}(undef, n)
		for i in 1:n
			symbols[i] = value_to_symbol_ordinal(series[i], breakpoints)
		end
		return symbols
	end
end

function discretize(delta::DeltaDiscretizer, series::Vector{Float64})
	diff_series = series
	for _ in 1:delta.order
		diff_series = diff(diff_series)
	end
	isnothing(delta.base_discretizer) || return discretize(delta.base_discretizer, diff_series)
	isempty(diff_series) && return Char[]
	min_val, max_val = extrema(diff_series)
	range_val = max_val - min_val
	if range_val < 1e-10
		middle_symbol = Char('a' + div(delta.alphabet_size, 2))
		return fill(middle_symbol, length(diff_series))
	end
	symbols = Vector{Char}(undef, length(diff_series))
	for i in 1:length(diff_series)
		normalized = (diff_series[i] - min_val) / range_val
		bin_idx = min(floor(Int, normalized * delta.alphabet_size), delta.alphabet_size - 1)
		symbols[i] = Char('a' + bin_idx)
	end
	return symbols
end

