


module TestMod

using Distributions
using MeasurementErrorModels.CovEst


function tacf(p::AbstractVector, lags::Int)
    a0, a1, a2 = p[1:3]
    ρ = zeros(lags)
    ρ[1] = 1.0
    ρ[2] = a1/(1.0-a2)
    for i in 3:(lags)
        ρ[i] = a1*ρ[i-1] + a2*ρ[i-2]
    end
    return ρ
end


function Asim(p::AbstractVector, lags::Int)
    ρ = tacf(p, lags)
    T11 = Matrix(CovEst.toeplitz(ρ, lags))
    T12 = zeros(lags, lags)
    ΣA = hvcat(2, T11, T12, T12, T11)
    g1 = MvNormal(zeros(2*lags), ΣA)
    return (A=reshape(rand(g1, 1), :, 2), ΣA=ΣA)
end


function Vsim(p::AbstractVector, lags::Int)
    ρ = tacf(p, lags)
    Σv = Matrix(CovEst.toeplitz(ρ, lags))
    g1 = MvNormal(zeros(lags), Σv)
    return (V=rand(g1, 1)[:], Σv=Σv)
end


end




