# MeasurementErrorModels.jl

[![Build Status](https://github.com/Mattriks/MeasurementErrorModels.jl/actions/workflows/CI.yml/badge.svg?branch=)](https://github.com/Mattriks/MeasurementErrorModels.jl/actions/workflows/CI.yml?query=branch%3A)
[![Coverage](https://codecov.io/gh/Mattriks/MeasurementErrorModels.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Mattriks/MeasurementErrorModels.jl)

## Summary

A julia package for working with Measurement Error Models.  This package was developed to calibrate Measurement Error Proxy System Models (MEPSMs).  
The package provides functions to initialise informative prior distributions, and a main function `Bmem()` to calculate a posterior distribution. The package uses [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) and [Measurements.jl](https://github.com/JuliaPhysics/Measurements.jl).

Inputs include:  
`df`: a dataframe with predictor variables and a response variable, given as [Measurement types](https://github.com/JuliaPhysics/Measurements.jl)  
`ΩV`: a 𝑛×𝑛 lag-covariance matrix for the response noise   
`ΣA`: a 𝑛𝑝×𝑛𝑝 lag-covariance matrix for the predictor noise   

An auxiliary function `dataf()` will provide defaults for `ΩV` and `ΣA` based on `df`.

```julia
using MeasurementErrorModels, LinearAlgebra, Distributions

μ₀, Σ₀ = [-0.22, 0.97*0.27], Diagonal([0.02^2, 0.15^2])
Bprior = Priors.makeprior(df, MvNormal(μ₀, Σ₀), ip=1000.0)

Y, ΩV, σV, X, ΣA, σA = dataf(df, ycol=ycol)
Bposterior = Bmem(Y, X, ΩV, ΣA, Bprior=Bprior)
```

Functions are also provided to estimate the noise auto- and cross-covariance in the predictors and response e.g.,

```julia
ΣA = CovEst.SigmaA(Bposterior.A, σA, lags=0:100)
```



# Installation

Install using the Julia [package manager](https://pkgdocs.julialang.org/v1/).  

From the Julia REPL, type `]` to enter the Pkg REPL mode and run:
```
pkg> add https://github.com/Mattriks/MeasurementErrorModels.jl
```
or from the Julia prompt:
```julia
julia> using Pkg; Pkg.add("https://github.com/Mattriks/MeasurementErrorModels.jl")
```



## Documentation

- [**LATEST**][docs-latest-url]


[docs-latest-url]: https://Mattriks.github.io/MeasurementErrorModels.jl/dev
[docs-stable-url]: https://Mattriks.github.io/MeasurementErrorModels.jl/stable
