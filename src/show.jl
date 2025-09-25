
function Base.show(io::IO, w::Word{D,T}) where {D,T}
	print(io, "Word{$(nameof(D)), $(T)}(")
 	print(io, String(w.symbols))
 	print(io, ")")
end


