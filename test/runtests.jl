using MeasurementErrorModels
using Test, Random
using CSV, DataFrames, Distributions
using LinearAlgebra
import Measurements as Mm
import ProfileLikelihood as pl
using OptimizationOptimJL

#@ @testset "MeasurementErrorModels.jl" begin
#    # Write your tests here.
# end


Munit =  Mm.Measurement{Float64}
f1 = joinpath(dirname(pathof(MeasurementErrorModels)), "../docs/src/assets/")
D2 = CSV.read(f1*"Palaui.csv", DataFrame, types=Dict([2,3,4].=>Munit))


@testset "Priors" begin
    B0 = MvNormal([-0.22, 0.27*0.97], [0.02^2 0.0; 0.0 0.15^2])
    Bprior = Priors.makeprior(D2[:,2:4], B0)
    @test mean(Bprior) ≈ [-8.38, -0.22, 0.2619] atol=0.01
end


@testset "CovEst" begin
    p = [0.0, 1.14, -0.22]
    lags = 200
    Random.seed!(1234)
    A, ΣA = TestMod.Asim(p, lags)
    σA = ones(lags, 2)
    SA = CovEst.SigmaA(A, σA, lags=0:50)
    @test norm(ΣA - SA, Inf)<0.5    

    Random.seed!(1234)
    V, Σv = TestMod.Vsim(p, lags)
    σV = ones(lags)
    Sv = CovEst.OmegaV(V, σV, lags=0:50)
    @test norm(Σv - Sv, Inf)<0.5
end


@testset "WTLS and HML" begin
    Bml = mean(wtls(D2[:,2:4]))
    data = dataf(D2[:,2:4], A="hannart")
    prob = pl.LikelihoodProblem(hloglik, [-10, -0.22, 0.27], data=data, syms=[:b₀, :b₁, :b₂])
    Bhml = pl.mle(prob, NelderMead()).mle
    @test Bml≈Bhml atol=0.001
end

