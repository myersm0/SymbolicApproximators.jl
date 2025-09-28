# SymbolicApproximators

## Overview
A Julia package implementing **Symbolic Aggregate approXimation (SAX)** and related methods for symbolic discretization and dimension reduction. These techniques allow you to take continuous-valued sequences (including but not necessarily limited to time series) and then use things like string algorithms and machine learning algorithms for categorical data.

The following algorithms are implemented or under development:
| |Algorithm|
|-|:-------------------------------------------------------------------------------------------------|
|☑|Piecewise Aggregate Approximation (PAA) ([Keogh et al 2000](https://www.cs.ucr.edu/~eamonn/SAX.pdf) and [Yi et al 2000](https://dl.acm.org/doi/10.5555/645926.671689))|
|☑|Symbolic Aggregate approXimation (SAX) ([Lin et al 2003](https://www.cs.ucr.edu/~eamonn/SAX.pdf))|
|☐|Extended SAX (ESAX) (Lkhagva et al 2006)|
|☐|Indexable SAX (iSAX) ([Shieh et al 2008](https://www.cs.ucr.edu/~eamonn/iSAX.pdf))|
|☐|Trend Feature SAX (TFSAX) ([Yu et al 2019](https://arxiv.org/abs/1905.00421)|
|☐|ordinal patterns ([Bandt & Pompe 2002](https://pubmed.ncbi.nlm.nih.gov/12005759/))|
|☐|delta encodings (todo: find reference)|

Note that the first one, PAA, is actually a continuous rather than a symbolic representation, but we include it in the package because all the SAX variants depend on it.

Coming soon there will be additional functionality such as:
- rolling window/segment functions to enable operations on streaming data
- distance functions (using [Distances.jl](https://github.com/JuliaStats/Distances.jl))
- other functions depending on algorithm, such as numerosity reduction and permutation entropy

## Usage
Usage revolves around this basic workflow:
1. Define a **`SymbolicApproximator`** with integer arguments for:
    - **word size**, i.e. the number of segments to use, or in other words the desired output length or dimensionality.
    - **alphabet size**, also called cardinality. This refers to the number of symbols (`Char`s by default) that you want to use.
        - Alternatively instead of specifying an alphabet _length_, you can directly pass in the _symbol set_ that you want to use in your output word. For example, `SAX(10, -2:2)` will give you an alphabet size of 5 where the symbols will be the integers -2 through 2 inclusive.
3. Pass that approximator and your data (presumably normalized -- see below) into the function **`encode(::SymbolicApproximator, ::AbstractVector)`**. Or if you prefer, use `approximate()` which is an alias for `encode()`.
4. Your output will be a **`Word`** composed of instances of the symbol set defined in the approximator.
    - A `Word` simply holds a vector of the symbols representing your output, along with a reference to the `SymbolicApproximator` that generated it (to inform things like the computation of distance between two `Word`s, for examle, which may be algorithm-dependent).
    - Since `Word` extends `AbstractVector`, you can simply use it as such.

Note that some algorithms expect preprocessed inputs. Specifically, SAX and variants expect the data to be normalized with a mean of 0, standard deviation of 1. It's left to the user to handle such preprocessing where necessary, mainly because there are a number of ways you can do it, depending on your situation: maybe your data already happens to be normally distributed in this manner, or maybe you have streaming data and need to do online normalization, etc.

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
