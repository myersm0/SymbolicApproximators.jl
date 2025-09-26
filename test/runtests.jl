using SymbolicApproximators
using StatsBase
using Test

@testset "SymbolicApproximators Tests" begin
	
	@testset "SAX basic functionality" begin
		# sine wave
		t = range(0, 2π, length=100)
		signal = sin.(t)
		
		# basic discretization
		sax = SAX(10, 3)  # 10 segments, 3 symbols (a,b,c)
		symbols = encode(sax, signal)
		@test length(symbols) == 10
		@test all(s in ['a', 'b', 'c'] for s in symbols)
		
		# test that sine wave produces expected pattern
		# (should go from low (a) to high (c) and back)
		@test symbols[1] in ['a', 'b']
		@test 'c' in symbols[3:7]  # has high values in middle
		@test symbols[end] in ['a', 'b']
		
		reconstructed = reconstruct(sax, symbols, 100)
		@test length(reconstructed) == 100
		# should preserve general shape of original
		@test cor(signal, reconstructed) > 0.7
	end
	
	@testset "SAX constant series handling" begin
		# Constant series should map to middle symbol
		constant_signal = ones(50)
		sax = SAX(5, 5)  # 5 symbols, so 'c' is middle
		symbols = discretize(sax, constant_signal)
		@test all(s == 'c' for s in symbols)
	end
	
	@testset "SAX distance measures" begin
		sax = SAX(8, 4)
		# identical sequences should have distance 0
		symbols1 = ['a', 'b', 'b', 'c', 'd', 'd', 'c', 'b']
		dist = distance(sax, symbols1, symbols1)
		@test dist ≈ 0.0
		# adjacent symbols should have distance 0
		symbols2 = ['b', 'c', 'c', 'd', 'd', 'c', 'b', 'a']
		dist = distance(sax, symbols1, symbols2)
		# all adjacent or same, so distance should be 0
		@test dist ≈ 0.0
		# very different sequences should have larger distance
		symbols3 = ['d', 'd', 'd', 'd', 'a', 'a', 'a', 'a']
		dist = distance(sax, symbols1, symbols3)
		@test dist > 0.0
	end
	
	@testset "OrdinalApproximator" begin
		# test ranking method
		data = [3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0]
		ord_rank = OrdinalApproximator(3, method=:rank)
		symbols = discretize(ord_rank, data)
		@test length(symbols) == length(data)
		@test all(s in ['a', 'b', 'c'] for s in symbols)
		# highest values should get highest symbols
		max_idx = argmax(data)
		@test symbols[max_idx] == 'c'
		# test quantile method
		ord_quant = OrdinalApproximator(3, method=:quantile)
		symbols = discretize(ord_quant, data)
		@test length(symbols) == length(data)
	end
	
	@testset "DeltaApproximator" begin
		linear = collect(1.0:10.0)
		# first-order differences should be constant
		delta = DeltaApproximator(3, order=1)
		symbols = discretize(delta, linear)
		@test length(symbols) == length(linear) - 1
		# all differences are the same (1.0), should map to same symbol
		@test all(s == symbols[1] for s in symbols)

		oscillating = [1.0, 3.0, 2.0, 4.0, 3.0, 5.0]
		symbols = discretize(delta, oscillating)
		# should alternate between positive and negative differences
		@test length(unique(symbols)) > 1
		# test with base discretizer
		delta_sax = DeltaApproximator(3, order=1, base=SAX(4, 3))
		t = range(0, 2π, length=50)
		signal = sin.(t)
		symbols = discretize(delta_sax, signal)
		@test length(symbols) == 4  # SAX reduces to 4 segments
	end
	
	@testset "Distance functions" begin
		# test that different discretizers have working distance functions
		data1 = collect(1:100) .+ 0.1
		data2 = collect(1:100) .* 2.0
		for disc in [SAX(10, 5), OrdinalApproximator(5), DeltaApproximator(5)]
			symbols1 = discretize(disc, data1)
			symbols2 = discretize(disc, data2)
			# distance to self should be 0
			@test distance(disc, symbols1, symbols1) ≈ 0.0
			# distance should be symmetric
			d12 = distance(disc, symbols1, symbols2)
			d21 = distance(disc, symbols2, symbols1)
			@test d12 ≈ d21
			# different sequences should have non-zero distance (usually)
			if symbols1 != symbols2
				@test distance(disc, symbols1, symbols2) > 0.0
			end
		end
	end
	
end

