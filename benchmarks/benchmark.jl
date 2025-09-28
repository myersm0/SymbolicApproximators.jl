
using SymbolicApproximators
using BenchmarkTools
using Chain

# simple test series with high-freq peaks
# (high-freq only on positive half)
signal = @chain begin
	range(0, 4π, length=500)
	sin.(_) + 0.5*sin.(20*_) .* (cos.(_) .> 0)  
	(_ .- mean(_)) ./ std(_)
end

@btime PAA(10)                    #     1 ns
model = PAA(10)
@btime encode($model, $signal)    #    86 ns
symbols = encode(model, signal)
@btime values($symbols)           #     1 ns
vals = values(symbols)

@btime SAX(10, 5)                 #   120 ns
model = SAX(10, 5)
@btime encode($model, $signal)    #   156 ns
symbols = encode(model, signal)
@btime values($symbols)           #   281 ns
vals = values(symbols)

@btime ESAX(10, 5)                #   118 ns
model = ESAX(10, 5)
@btime encode($model, $signal)    #  2782 ns
symbols = encode(model, signal)
@btime values($symbols)           # 24959 ns
vals = values(symbols)


signal = @chain begin
	range(0, 2π, length=500)
	sin.(_) + 0.5*sin.(20*_) .* (cos.(_) .> 0)  
	(_ .- mean(_)) ./ std(_)
end









