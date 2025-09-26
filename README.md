# SymbolicApproximators
A Julia package implementing Symbolic Aggregate approXimation (SAX) and related methods for symbolic time series discretization and dimension reduction.

In particular we have implemented (or are in the process of doing so) the following algorithms:
| |Algorithm|
|-|:------------------------------------------------------------------------------------------------|
|☑|Piecewise Aggregate Approximation (PAA) (todo: find reference)|
|☑|Symbolic Aggregate approXimation (SAX) ([Lin et al 2003](https://www.cs.ucr.edu/~eamonn/SAX.pdf))|
|☐|Extended SAX (ESAX) (Lkhagva et al 2006)|
|☐|Indexable SAX (iSAX) ([Shieh et al 2008](https://www.cs.ucr.edu/~eamonn/iSAX.pdf))|
|☐|Trend Feature SAX (TFSAX) ([Yu et al 2019](https://arxiv.org/abs/1905.00421)|
|☐|ordinal patterns ([Bandt & Pompe 2002](https://pubmed.ncbi.nlm.nih.gov/12005759/))|
|☐|delta encodings (todo: find reference)|

Note that the first one, PAA, is actually a continuous rather than a symbolic representation, but we include it in the package because all the SAX variants depend on it.

## Usage
```julia
using SymbolicApproximators
using StatsBase

signal = sin.(range(0, 4π, length=100))
normalized = (signal .- mean(signal)) ./ std(signal, corrected = false)

approximator = SAX(10, 5)  # discretize using SAX with 10 segments, 5 symbols
symbols = encode(approximator, normalized)
# Result: ['d', 'e', 'c', 'a', 'b', 'd', 'e', 'c', 'a', 'b']
```

## License

MIT


[![Build Status](https://github.com/myersm0/SymbolicApproximators.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/myersm0/SymbolicApproximators.jl/actions/workflows/CI.yml?query=branch%3Amain)
