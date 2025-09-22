
module SymbolicDiscretizers

using LinearAlgebra
using Distributions
using StatsBase

export SymbolicDiscretizer, SAX, OrdinalDiscretizer, DeltaDiscretizer
export discretize, reconstruct, distance
export StreamingDiscretizer, update!, get_symbols
export build_distance_table, numerosity_reduction

abstract type SymbolicDiscretizer end

struct SAX <: SymbolicDiscretizer
	nsegments::Int
	alphabet_size::Int
	normalize::Bool
	breakpoints::Vector{Float64}
	distance_table::Union{Nothing, Matrix{Float64}}
	function SAX(nsegments::Integer, alphabet_size::Integer; normalize = true, use_table = true)
		breakpoints = compute_gaussian_breakpoints(alphabet_size)
		distance_table = use_table ? build_distance_table(breakpoints) : nothing
		new(nsegments, alphabet_size, normalize, breakpoints, distance_table)
	end
end

struct OrdinalDiscretizer <: SymbolicDiscretizer
	alphabet_size::Int
	method::Symbol
	function OrdinalDiscretizer(alphabet_size::Int; method::Symbol=:quantile)
		method in [:rank, :quantile] || error("method must be :rank or :quantile")
		new(alphabet_size, method)
	end
end

struct DeltaDiscretizer <: SymbolicDiscretizer
	alphabet_size::Int
	order::Int
	base_discretizer::Union{Nothing, SymbolicDiscretizer}
	function DeltaDiscretizer(alphabet_size::Int; order::Int=1, base::Union{Nothing, SymbolicDiscretizer}=nothing)
		order > 0 || error("order must be positive")
		new(alphabet_size, order, base)
	end
end

mutable struct StreamingDiscretizer{T<:SymbolicDiscretizer}
	discretizer::T
	window_size::Int
	buffer::Vector{Float64}
	position::Int
	symbols::Vector{Char}
	function StreamingDiscretizer(disc::T, window_size::Int) where T<:SymbolicDiscretizer
		new{T}(disc, window_size, zeros(window_size), 0, Char[])
	end
end

function update!(sd::StreamingDiscretizer, value::Float64)
	sd.position = mod1(sd.position + 1, sd.window_size)
	sd.buffer[sd.position] = value
	# only discretize when we have a full window
	if sd.position == sd.window_size
		# create a properly ordered view of the circular buffer
		if sd.position == sd.window_size
			ordered_buffer = sd.buffer
		else
			ordered_buffer = vcat(sd.buffer[sd.position+1:end], sd.buffer[1:sd.position])
		end
		sd.symbols = discretize(sd.discretizer, ordered_buffer)
		return true
	end
	return false
end

function get_symbols(sd::StreamingDiscretizer)
	return sd.symbols
end

function build_distance_table(breakpoints::Vector{Float64})
	n = length(breakpoints) + 1
	table = zeros(n, n)
	for r in 1:n
		for c in 1:n
			if abs(r - c) <= 1
				table[r, c] = 0.0
			else
				max_idx = max(r, c)
				min_idx = min(r, c)
				if min_idx == 1
					low_val = breakpoints[1]
				else
					low_val = breakpoints[min_idx - 1]
				end
				if max_idx > length(breakpoints)
					high_val = breakpoints[end]
				else
					high_val = breakpoints[max_idx - 1]
				end
				table[r, c] = abs(high_val - low_val)
			end
		end
	end
	return table
end

function numerosity_reduction(symbols::Vector{Char})
	# run-length encoding for consecutive identical symbols
	vals, lens = rle(symbols)
	return vals, lens
end

function expand_numerosity(vals::Vector{Char}, lens::Vector{Int})
	return inverse_rle(vals, lens)
end


## helpers

function compute_gaussian_breakpoints(alphabet_size::Int)
	 alphabet_size > 1 || error("alphabet_size must be at least 2")
	 return quantile.(Normal(), (1:(alphabet_size-1)) ./ alphabet_size)
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

## distances

function distance(d::SAX, symbols1::Vector{Char}, symbols2::Vector{Char})
	length(symbols1) == length(symbols2) || error("Symbol sequences must have same length")
	if !isnothing(d.distance_table)
		# use precomputed lookup table
		dist_sum = 0.0
		for i in 1:length(symbols1)
			idx1 = Int(symbols1[i] - 'a') + 1
			idx2 = Int(symbols2[i] - 'a') + 1
			dist_sum += d.distance_table[idx1, idx2]^2
		end
		return sqrt(dist_sum)
	else
		# fall back to computation
		dist_sum = 0.0
		for i in 1:length(symbols1)
			dist_sum += distance(symbols1[i], symbols2[i], d.breakpoints)^2
		end
		return sqrt(dist_sum)
	end
end

function distance(s1::Char, s2::Char, breakpoints::Vector{Float64})
	s1 == s2 && return 0.0
	
	idx1 = Int(s1 - 'a')
	idx2 = Int(s2 - 'a')
	
	abs(idx1 - idx2) > 1 || return 0.0
	
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


## reconstruction

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

end
