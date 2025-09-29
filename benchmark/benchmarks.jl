
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

model = PAA(50)
symbols = encode(model, signal)
vals = values(symbols)
SUITE["PAA config"] = @benchmarkable PAA(10)
SUITE["PAA encode"] = @benchmarkable encode($model, $signal)
SUITE["PAA values"] = @benchmarkable values($symbols)

model = SAX(50, 25)
symbols = encode(model, signal)
vals = values(symbols)
SUITE["SAX config"] = @benchmarkable SAX(10, 5)
SUITE["SAX encode"] = @benchmarkable encode($model, $signal)
SUITE["SAX values"] = @benchmarkable values($symbols)

model = ESAX(50, 25)
symbols = encode(model, signal)
vals = values(symbols)
SUITE["ESAX config"] = @benchmarkable ESAX(10, 5)
SUITE["ESAX encode"] = @benchmarkable encode($model, $signal)
SUITE["ESAX values"] = @benchmarkable values($symbols)

run(SUITE)


