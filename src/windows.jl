
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

