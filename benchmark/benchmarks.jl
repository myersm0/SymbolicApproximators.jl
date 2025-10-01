
using SymbolicApproximators
using BenchmarkTools
using Chain

const SUITE = BenchmarkGroup()

# simple test series with high-freq peaks
# (high-freq only on positive half)
signal = @chain begin
	range(0, 4π, length = 1000)
	sin.(_) + 0.5*sin.(20*_) .* (cos.(_) .> 0)  
	(_ .- mean(_)) ./ std(_)
end


# 1 ns, 297 ns, 1 ns
model = PAA(50)
symbols = encode(model, signal)
symbols2 = encode(model, reverse(signal))
vals = values(symbols)
SUITE["PAA config"] = @benchmarkable PAA(50)
SUITE["PAA encode"] = @benchmarkable encode($model, $signal)
SUITE["PAA values"] = @benchmarkable values($symbols)
SUITE["PAA distance"] = @benchmarkable mindist($symbols, $symbols2)

# 225 ns, 688 ns, 24 ns
model = SAX(50, 25)
symbols = encode(model, signal)
symbols2 = encode(model, reverse(signal))
vals = values(symbols)
SUITE["SAX config"] = @benchmarkable SAX(50, 25)
SUITE["SAX encode"] = @benchmarkable encode($model, $signal)
SUITE["SAX values"] = @benchmarkable values($symbols)
SUITE["SAX distance"] = @benchmarkable mindist($symbols, $symbols2)

# 225 ns, 1906 ns, 94 ns
model = ESAX(50, 25)
symbols = encode(model, signal)
vals = values(symbols)
SUITE["ESAX config"] = @benchmarkable ESAX(50, 25)
SUITE["ESAX encode"] = @benchmarkable encode($model, $signal)
SUITE["ESAX values"] = @benchmarkable values($symbols)

run(SUITE)


