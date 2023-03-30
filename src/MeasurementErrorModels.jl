module MeasurementErrorModels

using LinearAlgebra, SparseArrays, Statistics
using DataFrames, Distributions
import Measurements as Mm

export Bmem, Bout, dataf, wtls, hloglik
export CovEst, Priors, TestMod






"""
    BBml(e, W::T Xa::T, ΣA::AbstractMatrix, BI::T) where T<:AbstractMatrix

Compute the covariance matrix of the weighted TLS solution (``\\Sigma_B``), using a Sandwich estimator.
Requires the computed error vector `e` and it's weight matrix `W`, the predictor matrix `Xa` (with intercept term) and its covariance matrix `ΣA`,
and ``BI = B ⊗ I_n``.
"""
function BBml(e, W::T, Xa::T, ΣA::AbstractMatrix, BI::T) where T<:AbstractMatrix
    n = length(e)
    Σ = Diagonal(e.^2)
    AW = U(e, W, ΣA, BI)
    m1 = (Xa'/W + AW)
    iXtX = inv(m1*Xa)
    M = m1 * Σ * m1'  # W or Σ?
    BB = (iXtX * (M) * (iXtX'))
    return 0.5*(BB+BB')
end


"""
    Bml(Y::AbstractVector, X::AbstractMatrix, Ωv::AbstractMatrix, ΣA::AbstractMatrix; Bprior=MvNormal(I(size(X,2)+1)))

Compute the weighted TLS solution {``B_\\text{ml}``}, using a fixed point iteration.  
Requires the response vector `Y` and the predictor matrix `X` (without intercept term), and their respective covariance matrices `Ωv` and `ΣA`.
"""
function Bml(Y::AbstractVector, X::AbstractMatrix, Ωv::AbstractMatrix, ΣA::AbstractMatrix; Bprior=MvNormal(I(size(X,2)+1)))
    n, p  = size(X)
    Xa = [ones(n) X]
    B = mean(Bprior)
    for i in 1:20
        BI = kron(B[2:end], I(n))
        W =  Ωv + BI'ΣA*BI
        e = W\(Y-Xa*B)
        G = kron(I(p), e')*ΣA*BI/W
        U = sparse([zeros(1, n); G])
        B = (Xa'/W*Xa + U*Xa) \ ((Xa'/W + U)*Y)
#        @show B
    end
    return B
end


"""
    Bmem(Y::AbstractVector, X::AbstractMatrix, Ωv::AbstractMatrix, ΣA::AbstractMatrix; Bprior=MvNormal(I(size(X,2)+1)))

Compute a quasi-Bayesian posterior for the Measurement Error Model (one iteration).
Requires the response vector `Y` and the predictor matrix `X` (without intercept term), and their respective covariance matrices `Ωv` and `ΣA`.
"""
function Bmem(Y::AbstractVector, X::AbstractMatrix, Ωv::AbstractMatrix, ΣA::AbstractMatrix; Bprior=MvNormal(I(size(X,2)+1)))
    n, p  = size(X)
    Xa = [ones(n) X]
    Bpr, QB = mean(Bprior), cov(Bprior)
    PB = inv(QB)
# Maximum likelihood
    Bml1 = Bml(Y, X, Ωv, ΣA)
    e = Y - Xa*Bml1
    BI = kron(Bml1[2:end], I(n))
    W =  Ωv + BI'ΣA*BI
    BBmli = inv(BBml(e, W, Xa, ΣA, BI))
# Bayesian posterior
    B = (BBmli + PB)\(BBmli*Bml1 + PB*Bpr)
    covB = inv(BBmli + PB)
    e = Y - Xa*B
    BI = kron(B[2:end], I(n))
    W =  Ωv + BI'ΣA*BI
    A = Af(e, W, ΣA, BI)
    V = Vf(e, W, Ωv)
    return (e=e, V=V, A=A, Bpost=MvNormal(B,  0.5*(covB+covB')))
end




function Af(e, W::AbstractMatrix, ΣA::AbstractMatrix, BI::AbstractMatrix)
    p = size(BI,1) ÷ size(BI,2)
    A = -kron(I(p), (W\e)')*ΣA*BI
    return permutedims(A)
end


function Vf(e, W::AbstractMatrix, ΩV::AbstractMatrix)
    V = ΩV*(W\e)
    return V
end


function U(e, W::AbstractMatrix, ΣA::AbstractMatrix, BI::AbstractMatrix)
    # U = AW
    n = size(BI,2)
    p = size(BI,1) ÷ n
    G = kron(I(p), (W\e)')*ΣA*BI/W
    AW = sparse([zeros(1, n); G])
    return AW
end


"""
    Bout(D::DataFrame; Bprior=MvNormal(size(D,2)-1), lags=0:25, ycol=:ya)

Compute a quasi-Bayesian posterior for the Measurement Error Model (two-step iteration).
Requires a DataFrame `D` containing all variables, the name of the response variable can be specified by `ycol`.
`lags` are the lags retained for modelling the lagged covariance matrix of the predictor noise, ``\\Sigma_A``.
"""
function Bout(D::DataFrame; Bprior=MvNormal(I(2)), lags=0:25, ycol=:ya)
    Y, ΩV, σV, X, ΣA, σA = dataf(D, ycol=ycol)
    Bw1 = Bmem(Y, X, ΩV, ΣA, Bprior=Bprior)
    ΣA = CovEst.SigmaA(Bw1.A, σA, lags=lags)
    Bw2 = Bmem(Y, X, ΩV, ΣA, Bprior=Bprior)
    return DataFrame(Dist=["Bprior", "Bpost1", "Bpost2"] , value=[Bprior, Bw1.Bpost, Bw2.Bpost])
end


function dataf(D::DataFrame; ycol=:ya, A="full")
    Y0, X0 = D[:,ycol], Matrix(D[:, Not(ycol)])
    X, Y = Mm.value.(X0), Mm.value.(Y0)
    σA, σV = Mm.uncertainty.(X0), Mm.uncertainty.(Y0)
    varA, varV = σA.^2, σV.^2
    ΣA = (A=="full") ? spdiagm(vec(varA)) : [spdiagm(vA) for vA in eachcol(varA)]
    ΩV = spdiagm(varV)
    return (Y=Y, ΩV=ΩV, σV=σV, X=X, ΣA=ΣA, σA=σA)
end


function wtls(D::DataFrame; ycol=:ya)
    Y, Ωv, σV, X, ΣA, σA = dataf(D, ycol=ycol)    
    n, p  = size(X)
    Xa = [ones(n) X]
#   Maximum likelihood
    Bml1 = Bml(Y, X, Ωv, ΣA)
    e = Y - Xa*Bml1
    BI = kron(Bml1[2:end], I(n))
    W =  Ωv + BI'ΣA*BI
    covB = BBml(e, W, Xa, ΣA, BI)
    return MvNormal(Bml1,  covB)
end



include("priors.jl")
include("covest.jl")
include("testmod.jl")


########## Maximum-likehood solution of Hannart et al. (2014) (modified)  ##########


function hloglik(θ, data)
    Y, Ωv, _, X, Ωa, _ = data
    p = size(X,2)
    b₀, B = θ[1], θ[2:end]
    Σa = sum([B[i]*B[i]*Ωa[i] for i in 1:p])
    e = Y .- b₀ .- X*B
    We = (Ωv + Σa)\e
    A = -reduce(hcat, [B[i]*Ωa[i]*We for i in 1:p])
    Xstar = X - A
    e1 = Y .-b₀ .- Xstar*B
    ll = e1'/Ωv*e1 + sum([A[:,i]'/Ωa[i]*A[:,i] for i in 1:p])
    return -0.5*ll
end
















end
