using MeasurementErrorModels
using Test
using CSV, DataFrames, Distributions
import Measurements as Mm

#@ @testset "MeasurementErrorModels.jl" begin
#    # Write your tests here.
# end


@testset "Priors" begin
    Munit =  Mm.Measurement{Float64}
    f1 = joinpath(dirname(pathof(MeasurementErrorModels)), "../docs/src/assets/")
    D2 = CSV.read(f1*"Palaui.csv", DataFrame, types=Dict([2,3,4].=>Munit))

    B0 = MvNormal([-0.22, 0.27*0.97], [0.02^2 0.0; 0.0 0.15^2])
    Bprior = Priors.makeprior(D2[:,2:4], B0)
    @test mean(Bprior) ≈ [-8.38, -0.22, 0.2619] atol=0.01
end



