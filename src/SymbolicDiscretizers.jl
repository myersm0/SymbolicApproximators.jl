
module SymbolicDiscretizers

using Statistics: mean, std, quantile
using LinearAlgebra: norm

export SymbolicDiscretizer, SAX, OrdinalDiscretizer, DeltaDiscretizer
export discretize, reconstruct, distance

abstract type SymbolicDiscretizer end

struct SAX <: SymbolicDiscretizer
	nsegments::Int  # dimensionality reduction (w in paper)
	alphabet_size::Int
	normalize::Bool
	breakpoints::Vector{Float64}
	function SAX(nsegments::Int, alphabet_size::Int; normalize::Bool=true)
		breakpoints = compute_gaussian_breakpoints(alphabet_size)
		new(nsegments, alphabet_size, normalize, breakpoints)
	end
end

struct OrdinalDiscretizer <: SymbolicDiscretizer
	alphabet_size::Int
	method::Symbol  # :rank or :quantile for now
	function OrdinalDiscretizer(alphabet_size::Int; method::Symbol=:quantile)
		method in [:rank, :quantile] || error("method must be :rank or :quantile")
		new(alphabet_size, method)
	end
end

struct DeltaDiscretizer <: SymbolicDiscretizer
	alphabet_size::Int
	order::Int  # 1 for first difference, 2 for second, etc.
	base_discretizer::Union{Nothing, SymbolicDiscretizer}
	function DeltaDiscretizer(alphabet_size::Int; order::Int=1, base::Union{Nothing, SymbolicDiscretizer}=nothing)
		order > 0 || error("order must be positive")
		new(alphabet_size, order, base)
	end
end


## helpers

function compute_gaussian_breakpoints(alphabet_size::Int)
	alphabet_size > 1 || error("alphabet_size must be at least 2")
	breakpoint_table = Dict(
		2 => Float64[0.0],
		3 => [-0.43, 0.43],
		4 => [-0.67, 0.0, 0.67],
		5 => [-0.84, -0.25, 0.25, 0.84],
		6 => [-0.97, -0.43, 0.0, 0.43, 0.97],
		7 => [-1.07, -0.57, -0.18, 0.18, 0.57, 1.07],
		8 => [-1.15, -0.67, -0.32, 0.0, 0.32, 0.67, 1.15],
	)
	haskey(breakpoint_table, alphabet_size) && return breakpoint_table[alphabet_size]
	# compute breakpoints for equal-area under standard normal
	breakpoints = Float64[]
	for i in 1:(alphabet_size-1)
		prob = i / alphabet_size
		push!(breakpoints, quantile_normal(prob))
	end
	return breakpoints
end

# todo: do this better with Distributions.jl
function quantile_normal(p::Float64)
	if p == 0.5
		return 0.0
	elseif p < 0.5
		return -sqrt(2) * erfcinv(2 * p)
	else
		return sqrt(2) * erfcinv(2 * (1 - p))
	end
end

# approximation of erfcinv for self-contained implementation
erfcinv(x::Float64) = -invnorm((x - 1) / 2) / sqrt(2)
invnorm(p::Float64) = sqrt(2) * erfinv(2 * p - 1)

# basic approximation of erfinv
function erfinv(x::Float64)
	a = 0.147
	sgn = sign(x)
	x = abs(x)
	lnx = log(1 - x^2)
	tt1 = 2 / (π * a) + lnx / 2
	tt2 = lnx / a
	return sgn * sqrt(sqrt(tt1^2 - tt2) - tt1)
end



# PAA (Piecewise Aggregate Approximation)
function paa(series::Vector{Float64}, nsegments::Int)
	n = length(series)
	segment_length = n / nsegments
	result = zeros(nsegments)
	for i in 1:nsegments
		start_idx = round(Int, (i - 1) * segment_length + 1)
		end_idx = round(Int, i * segment_length)
		end_idx = min(end_idx, n)
		result[i] = mean(series[start_idx:end_idx])
	end
	return result
end


## discretization

function discretize(sax::SAX, series::Vector{Float64})
	n = length(series)
	
	if sax.normalize
		μ = mean(series)
		σ = std(series, corrected=false)
		# Handle constant series
		if σ < 1e-10
			# Return middle symbol for constant series
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

function value_to_symbol(value::Float64, breakpoints::Vector{Float64})
	for i in 1:length(breakpoints)
		if value < breakpoints[i]
			return Char('a' + i - 1)
		end
	end
	return Char('a' + length(breakpoints))
end

# todo: use traits
function discretize(ord::OrdinalDiscretizer, series::Vector{Float64})
	n = length(series)
	if ord.method == :rank
		ranks = sortperm(sortperm(series))  # double sortperm gives ranks
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

function value_to_symbol_ordinal(value::Float64, breakpoints::Vector{Float64})
	for i in 1:length(breakpoints)
		if value <= breakpoints[i]
			return Char('a' + i - 1)
		end
	end
	return Char('a' + length(breakpoints))
end

function linspace_bins(start, stop, n)
	return collect(range(start, stop, length=n+1))[2:end]
end

function discretize(delta::DeltaDiscretizer, series::Vector{Float64})
	diff_series = series
	for _ in 1:delta.order
		diff_series = diff(diff_series)
	end
	if delta.base_discretizer !== nothing
		# Use base discretizer on the differences
		return discretize(delta.base_discretizer, diff_series)
	else
		# Simple equal-width binning of differences
		if isempty(diff_series)
			return Char[]
		end
		min_val, max_val = extrema(diff_series)
		range_val = max_val - min_val
		if range_val < 1e-10
			# All differences are the same
			middle_symbol = Char('a' + div(delta.alphabet_size, 2))
			return fill(middle_symbol, length(diff_series))
		end
		symbols = Vector{Char}(undef, length(diff_series))
		for i in 1:length(diff_series)
			# Map to [0, 1], then to alphabet
			normalized = (diff_series[i] - min_val) / range_val
			bin_idx = min(floor(Int, normalized * delta.alphabet_size), delta.alphabet_size - 1)
			symbols[i] = Char('a' + bin_idx)
		end
		return symbols
	end
end


## distances

function distance(d::SAX, symbols1::Vector{Char}, symbols2::Vector{Char})
	length(symbols1) == length(symbols2) || error("Symbol sequences must have same length")
	dist_sum = 0.0
	for i in 1:length(symbols1)
		dist_sum += distance(symbols1[i], symbols2[i], d.breakpoints)^2
	end
	return sqrt(dist_sum) # scale by compression ratio (as in SAX paper MINDIST)
end

function distance(s1::Char, s2::Char, breakpoints::Vector{Float64})
	s1 == s2 && return 0.0
	
	idx1 = Int(s1 - 'a')
	idx2 = Int(s2 - 'a')
	
	abs(idx1 - idx2) > 1 || return 0.0
	
	# get the breakpoint values
	max_idx = max(idx1, idx2)
	min_idx = min(idx1, idx2)
	
	if min_idx == 0
		low_val = breakpoints[1]
	else
		low_val = breakpoints[min_idx]
	end
	
	if max_idx > length(breakpoints)
		high_val = breakpoints[end]
	else
		high_val = breakpoints[max_idx]
	end
	
	return abs(high_val - low_val)
end

function distance(d::OrdinalDiscretizer, symbols1::Vector{Char}, symbols2::Vector{Char})
	# simple Hamming distance for ordinal
	length(symbols1) == length(symbols2) || error("Symbol sequences must have same length")
	return sum(symbols1 .!= symbols2)
end

function distance(d::DeltaDiscretizer, symbols1::Vector{Char}, symbols2::Vector{Char})
	# Hamming distance for delta (could be refined)
	length(symbols1) == length(symbols2) || error("Symbol sequences must have same length")
	return sum(symbols1 .!= symbols2)
end


## reconstruction (approximate inverse)

function reconstruct(sax::SAX, symbols::Vector{Char}, original_length::Int)
	# map symbols back to breakpoint midpoints
	paa_values = zeros(length(symbols))
	
	for i in 1:length(symbols)
		idx = Int(symbols[i] - 'a')
		if idx == 0
			# below first breakpoint
			paa_values[i] = sax.breakpoints[1] - 0.5
		elseif idx >= length(sax.breakpoints)
			# above last breakpoint  
			paa_values[i] = sax.breakpoints[end] + 0.5
		else
			# between breakpoints
			paa_values[i] = (sax.breakpoints[idx] + sax.breakpoints[idx+1]) / 2
		end
	end
	
	# inverse PAA - simple repetition
	segment_length = original_length / length(symbols)
	reconstructed = Float64[]
	
	for i in 1:length(symbols)
		repeat_count = round(Int, segment_length)
		if i == length(symbols)
			# last segment gets remaining points
			repeat_count = original_length - length(reconstructed)
		end
		append!(reconstructed, fill(paa_values[i], repeat_count))
	end
	
	return reconstructed[1:original_length]
end

end

