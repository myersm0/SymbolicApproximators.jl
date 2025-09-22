
abstract type SymbolicDiscretizer end

struct SAX <: SymbolicDiscretizer
	nsegments::Int
	alphabet_size::Int
	normalize::Bool
	breakpoints::Vector{Float64}
	distance_table::Union{Nothing, Matrix{Float64}}
	function SAX(nsegments::Integer, alphabet_size::Integer; normalize = true, use_table = true)
		alphabet_size > 1 || error("alphabet_size must be at least 2")
		function compute_gaussian_breakpoints(alphabet_size)
			return quantile.(Normal(), (1:(alphabet_size-1)) ./ alphabet_size)
		end
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

