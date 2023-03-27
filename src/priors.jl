

module Priors


using LinearAlgebra, Statistics
using DataFrames, Distributions
import Measurements as Mm




"""
    makeprior(Y::AbstractVector, X::AbstractMatrix, Bprior=MvNormal(I(size(D,2)-1)); ip=1.0)

Construct a prior with an informed intercept from `Y`, `X` and `Bprior`, 
where `Bprior` is the distribution of the predictor coefficients (not including the intercept).
A prior for the intercept is calculated using  `Y`, `X` and `Bprior`.
`ip` scales the estimated variance of the intercept.
"""
function makeprior(Y::AbstractVector, X::AbstractMatrix, Bprior=MvNormal(I(size(D,2)-1)); ip=1.0)
    n = length(Y)
    Bpr = mean(Bprior)
    μX = mean(X, dims=1)
    b₀ = mean(Y) .- μX*Bpr
    B = [b₀; Bpr]
    yh = [ones(n) X]*B
    σ² = var(Y-yh)
    RSS₀ = n .- (n^2)*(μX/ (X'X) * μX') 
    print("Var(b₀) = ", round.([σ² RSS₀], digits=4), "\n")
    QB = cat(ip*σ²/RSS₀[1], cov(Bprior), dims=(1,2))
#    PB = inv(QB)
    return MvNormal(B, QB)
end



"""
    makeprior(D::DataFrame, Bprior=MvNormal(I(size(D,2)-1)); ycol=:ya, ip=1.0)

Construct a prior with an informed intercept from a dataframe `D`. `ycol` is the column name of the response variable.
`Bprior` is the distribution of the predictor coefficients (not including the intercept).
`ip` scales the estimated variance of the intercept.
"""
function makeprior(D::DataFrame, Bprior=MvNormal(I(size(D,2)-1)); ycol=:ya, ip=1.0)
  X, Y = Matrix(Mm.value.(D[:, Not(ycol)])), Mm.value.(D[:,ycol])
  z =  makeprior(Y, X, Bprior, ip=ip)
  return z
end


end


