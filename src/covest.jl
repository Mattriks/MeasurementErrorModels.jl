




module CovEst

using LinearAlgebra, ToeplitzMatrices
using StatsBase: autocor, crosscor



function toeplitz(r, n::Int)
    nr = length(r)
    v = [r; zeros(n-nr)]
    M = Toeplitz(v,v)
    return M
end

"""
    OmegaV(V::AbstractVector, εV::AbstractVector; lags=0:25)

Estimate covariance matrix ``\\Omega_V`` from the error vector `V`, subject to known `εV`.
``n``-length vector `εV` is the heteroscedastic standard deviation of ``V``.
The underlying matrix Ω is assumed to be toeplitz for lags `lags`.
"""
function OmegaV(V::AbstractVector, εV::AbstractVector; lags=0:25)
    n = size(V,1)
    acf = autocor(V, lags)
    R = Matrix(toeplitz(acf, n))
    Rm = nearPD(0.5*(R+transpose(R)))
    Σv = Diagonal(εV)
    Ωy = Σv'*Rm*Σv
    return 0.5*(Ωy+Ωy')
end


"""
    SigmaA(A::AbstractMatrix, εA::AbstractMatrix; lags=0:25)

Estimate ``\\Sigma_A`` from the error field `A`, subject to known `εA`.
The underlying Ω matrices are assumed to be toeplitz for lags `lags`. 

"""
function SigmaA(A::AbstractMatrix, εA::AbstractMatrix; lags=0:25)
    n, p = size(A)
    ccf = crosscor(A, A, lags)
    Ωs = [Matrix(toeplitz(ccf[:,j,k], n))  for j in axes(ccf,2), k in axes(ccf,3)] 
    # R = cross-correlation matrix
    R = hvcat(p, Ωs...)
    Rm = nearPD(0.5*(R+transpose(R)))
    Σa = Diagonal(vec(εA))
    ΣA = Σa'*Rm*Σa
    return 0.5*(ΣA+ΣA')
end



function nearPD(M::AbstractMatrix; ep=eps(Float64))
    L, V = eigen(M)
    BB = V*Diagonal(max.(L, 10*ep))*(V')
    T = 1.0./sqrt.(diag(BB))
    TT = T .* T'
    z = BB .* TT
    return 0.5*(z+z') 
end





end














