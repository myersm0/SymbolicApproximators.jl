
function Base.show(io::IO, ::MIME"text/plain", ca::ContinuousApproximator)
	name = typeof(ca).name.name
	w = word_size(ca)
	print(io, name, "(", w, ")")
end

function Base.show(io::IO, ::MIME"text/plain", sa::SymbolicApproximator)
	name = typeof(sa).name.name
	w = word_size(sa)
	c = cardinality(sa)
	print(io, name, "(", w, ", ", c, ")")
end

function Base.show(io::IO, ::MIME"text/plain", w::Word{<:SymbolicApproximator})
	approx = w.approximator[]
	name = typeof(approx).name.name
	symbols = values(w)
	print(io, name, " Word: ")
	if eltype(symbols) <: Char
		print(io, "\"", join(symbols), "\"")
	else
		print(io, symbols)
	end
end

function Base.show(io::IO, ::MIME"text/plain", w::Word{<:ContinuousApproximator})
	approx = w.approximator[]
	name = typeof(approx).name.name
	vals = values(w)
	print(io, name, " Word: ", round.(vals, digits=3))
end

