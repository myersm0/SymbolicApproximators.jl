
using SymbolicApproximators
using StatsBase
using DataStructures
using Observables
using GLMakie

w = 10             # word size
a = 5              # alphabet size
window_len = 200   # length of time series subsequeences

# visualization refresh rate
fps = 30
Δt = 1 / fps

sax = SAX(w, a)

# allocate a buffer for input accumulation
buf = CircularBuffer{Float64}(window_len)

# pre-allocate a vector to store SAX words,
# to avoid allocating new memory at every iteration
dest = Vector{Int}(undef, w)

# observables for reactive updates
signal = Observable(Vector{Float64}(undef, window_len))
words = Observable(Vector{Int}(undef, w))

# visualization setup
fig = Figure(resolution = (800, 400))
ax1 = Axis(fig[1, 1], title = "Incoming signal")
lineplot = lines!(ax1, 1:window_len, signal)
ax2 = Axis(
	fig[2, 1], 
	title = "Current SAX word (symbol indices)", 
	limits = (1, w, 0, a + 1)
)

barplot = barplot!(ax2, 1:w, words)
display(fig)

@async begin
	t = 0.0
	while true
		# simulate new sample
		new_val = sin(t) + 0.1randn() # noisy sine wave
		push!(buf, new_val)
		if length(buf) == window_len
			# normalize
			normed = (buf .- mean(buf)) ./ std(buf, corrected=false)
			# encode into preallocated vector (mutating)
			w = encode!(sax, dest, normed)
			# update observables
			signal[] = copy(buf)		  # make a safe copy for plotting
			words[] = copy(keys(w))	  # make a safe copy (ephemeral)
		end
		sleep(1 / fps)
		t += 0.1
	end
end




