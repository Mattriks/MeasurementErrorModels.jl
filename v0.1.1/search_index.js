var documenterSearchIndex = {"docs":
[{"location":"man/Example2/#Example-2","page":"Example 2","title":"Example 2","text":"","category":"section"},{"location":"man/Example2/","page":"Example 2","title":"Example 2","text":"In this example, the first posterior (Bpost1 in the table below) \"approaches\" the basic ML solution (see Example 3), except the prior has λ=10. The final posterior (Bpost2) includes the full auto- and cross-correlation in Σ_A, so is different from the basic ML solution.","category":"page"},{"location":"man/Example2/","page":"Example 2","title":"Example 2","text":"using CSV, DataFrames, Distributions, LinearAlgebra\nimport Measurements as Mm\nusing Measurements: ±\nusing MeasurementErrorModels\n\nMunit =  Mm.Measurement{Float64}\nf1 = joinpath(dirname(pathof(MeasurementErrorModels)), \"../docs/src/assets/\")\nD2 = CSV.read(f1*\"Palaui.csv\", DataFrame, types=Dict([2,3,4].=>Munit))\nfirst(D2, 3)","category":"page"},{"location":"man/Example2/","page":"Example 2","title":"Example 2","text":"\nfunction marginals(D::AbstractVector)\n    z = mean(D[1]) .± sqrt.(var(D[1]))\n    return (b₀=z[1], b₁=z[2], b₂=z[3])\nend\n\nμ₀, Σ₀ = [-0.22, 0.97*0.27], Diagonal([0.02^2, 0.15^2])\nBprior = Priors.makeprior(D2[:,2:4], MvNormal(μ₀, Σ₀), ip=1000.0)\n\nDbs = Bout(D2[:,2:end], Bprior=Bprior, lags=0:100)\nD3 = combine(groupby(Dbs, :Dist), :value=>marginals=>AsTable)\nD3","category":"page"},{"location":"lib/library/#Library","page":"Library","title":"Library","text":"","category":"section"},{"location":"lib/library/#Contents","page":"Library","title":"Contents","text":"","category":"section"},{"location":"lib/library/","page":"Library","title":"Library","text":"Modules = [MeasurementErrorModels, Priors, CovEst]","category":"page"},{"location":"lib/library/#Functions","page":"Library","title":"Functions","text":"","category":"section"},{"location":"lib/library/","page":"Library","title":"Library","text":"Modules = [MeasurementErrorModels, Priors, CovEst]","category":"page"},{"location":"lib/library/#MeasurementErrorModels.BBml-Union{Tuple{T}, Tuple{Any, T, T, AbstractMatrix, T}} where T<:(AbstractMatrix)","page":"Library","title":"MeasurementErrorModels.BBml","text":"BBml(e, W::T Xa::T, ΣA::AbstractMatrix, BI::T) where T<:AbstractMatrix\n\nCompute the covariance matrix of the weighted TLS solution (Sigma_B), using a Sandwich estimator. Requires the computed error vector e and it's weight matrix W, the predictor matrix Xa (with intercept term) and its covariance matrix ΣA, and BI = B  I_n.\n\n\n\n\n\n","category":"method"},{"location":"lib/library/#MeasurementErrorModels.Bmem-Tuple{AbstractVector, AbstractMatrix, AbstractMatrix, AbstractMatrix}","page":"Library","title":"MeasurementErrorModels.Bmem","text":"Bmem(Y::AbstractVector, X::AbstractMatrix, Ωv::AbstractMatrix, ΣA::AbstractMatrix; Bprior=MvNormal(I(size(X,2)+1)))\n\nCompute a quasi-Bayesian posterior for the Measurement Error Model (one iteration). Requires the response vector Y and the predictor matrix X (without intercept term), and their respective covariance matrices Ωv and ΣA.\n\n\n\n\n\n","category":"method"},{"location":"lib/library/#MeasurementErrorModels.Bml-Tuple{AbstractVector, AbstractMatrix, AbstractMatrix, AbstractMatrix}","page":"Library","title":"MeasurementErrorModels.Bml","text":"Bml(Y::AbstractVector, X::AbstractMatrix, Ωv::AbstractMatrix, ΣA::AbstractMatrix; Bprior=MvNormal(I(size(X,2)+1)))\n\nCompute the weighted TLS solution {B_textml}, using a fixed point iteration.   Requires the response vector Y and the predictor matrix X (without intercept term), and their respective covariance matrices Ωv and ΣA.\n\n\n\n\n\n","category":"method"},{"location":"lib/library/#MeasurementErrorModels.Bout-Tuple{DataFrames.DataFrame}","page":"Library","title":"MeasurementErrorModels.Bout","text":"Bout(D::DataFrame; Bprior=MvNormal(size(D,2)-1), lags=0:25, ycol=:ya)\n\nCompute a quasi-Bayesian posterior for the Measurement Error Model (two-step iteration). Requires a DataFrame D containing all variables, the name of the response variable can be specified by ycol. lags are the lags retained for modelling the lagged covariance matrix of the predictor noise, Sigma_A.\n\n\n\n\n\n","category":"method"},{"location":"lib/library/#MeasurementErrorModels.Priors.makeprior","page":"Library","title":"MeasurementErrorModels.Priors.makeprior","text":"makeprior(Y::AbstractVector, X::AbstractMatrix, Bprior=MvNormal(I(size(D,2)-1)); ip=1.0)\n\nConstruct a prior with an informed intercept from Y, X and Bprior,  where Bprior is the distribution of the predictor coefficients (not including the intercept). A prior for the intercept is calculated using  Y, X and Bprior. ip scales the estimated variance of the intercept.\n\n\n\n\n\n","category":"function"},{"location":"lib/library/#MeasurementErrorModels.Priors.makeprior-2","page":"Library","title":"MeasurementErrorModels.Priors.makeprior","text":"makeprior(D::DataFrame, Bprior=MvNormal(I(size(D,2)-1)); ycol=:ya, ip=1.0)\n\nConstruct a prior with an informed intercept from a dataframe D. ycol is the column name of the response variable. Bprior is the distribution of the predictor coefficients (not including the intercept). ip scales the estimated variance of the intercept.\n\n\n\n\n\n","category":"function"},{"location":"lib/library/#MeasurementErrorModels.CovEst.Omegav-Tuple{AbstractVector, AbstractVector}","page":"Library","title":"MeasurementErrorModels.CovEst.Omegav","text":"Omegav(A::AbstractMatrix, εV::AbstractMatrix; lags=0:25)\n\nEstimate covariance matrix Omega_v from the error vector V, subject to known εV. n-length vector εV is the heteroscedastic standard deviation of V. The underlying matrix Ω is assumed to be toeplitz for lags lags.\n\n\n\n\n\n","category":"method"},{"location":"lib/library/#MeasurementErrorModels.CovEst.SigmaA-Tuple{AbstractMatrix, AbstractMatrix}","page":"Library","title":"MeasurementErrorModels.CovEst.SigmaA","text":"SigmaA(A::AbstractMatrix, εA::AbstractMatrix; lags=0:25)\n\nEstimate Sigma_A from the error field A, subject to known εA. The underlying Ω matrices are assumed to be toeplitz for lags lags. \n\n\n\n\n\n","category":"method"},{"location":"theory/#Theory","page":"Theory","title":"Theory","text":"","category":"section"},{"location":"theory/#General-EIV-model","page":"Theory","title":"General EIV model","text":"","category":"section"},{"location":"theory/","page":"Theory","title":"Theory","text":"beginequation\nbeginaligned\n    vb*z_t= vb*z^*_t + vb*q_t \n    vb*z^*_t = left y^*_t quad vb*x^*_tright \n    vb*q_t = leftv_t quad  vb*a_t right \n    y_t^* = b_0 + vb*x^*_tvb*B^* + epsilon_t\nendaligned\nendequation","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"beginalign\nbeginbmatrix Z_11  Z_21  vdots   Z_n1endbmatrix = \n    beginbmatrix\n        1  vb* x_1  1  vb* x_2  vdots  1  vb* x_n \n    endbmatrix\n    beginbmatrix b_0  vb*B^*    endbmatrix + \n    beginbmatrixmathbbI_n quad -(vb*B^* otimes mathbbI_n)    endbmatrix\n    beginbmatrix vb*V  vec(mathbf A)   endbmatrix\nendalign","category":"page"},{"location":"theory/#Maximum-likelihood","page":"Theory","title":"Maximum likelihood","text":"","category":"section"},{"location":"theory/","page":"Theory","title":"Theory","text":"beginalign\nvb B_textwtls =left (vb X-vb A)vb W^-1vb X right ^-1left(vb X-vb A)vb W^-1vb Yright    \nSigma_B^textwtls =left (vb X-vb A)vb W^-1vb X right^-1mathbfMleft (vb X-vb A)vb W^-1vb X right^-1\nendalign","category":"page"},{"location":"theory/","page":"Theory","title":"Theory","text":"where mathbfM=(vb X-vb A)vb W^-1Sigma vb W^-1(vb X-vb A).","category":"page"},{"location":"theory/#Bayesian-solution","page":"Theory","title":"Bayesian solution","text":"","category":"section"},{"location":"theory/","page":"Theory","title":"Theory","text":"beginalign\n    vb B = left(Sigma_B^textwtls^-1 + Sigma_0^-1right)^-1left(Sigma_B^textwtls^-1vb B_textwtls+Sigma_0^-1vb B_0right) \n    Sigma_B =left( Sigma_B^textwtls^-1+Sigma_0^-1right)^-1\nendalign","category":"page"},{"location":"theory/#Reference","page":"Theory","title":"Reference","text":"","category":"section"},{"location":"theory/","page":"Theory","title":"Theory","text":"Fischer, M (2023) Measurement Error Proxy System Models: MEPSM v0.2. link to preprint","category":"page"},{"location":"man/Example1/#Example1","page":"Example 1","title":"Example 1","text":"","category":"section"},{"location":"man/Example1/","page":"Example 1","title":"Example 1","text":"In this example, the final posterior is calculated for a range of λ values.   The solution for λ=10 is the same as Bpost2 in Example 2.","category":"page"},{"location":"man/Example1/","page":"Example 1","title":"Example 1","text":"\nusing CSV, DataFrames, Distributions, LinearAlgebra\nimport Measurements as Mm\nusing MeasurementErrorModels\n\nMunit =  Mm.Measurement{Float64}\nf1 = joinpath(dirname(pathof(MeasurementErrorModels)), \"../docs/src/assets/\")\nD2 = CSV.read(f1*\"Palaui.csv\", DataFrame, types=Dict([2,3,4].=>Munit))\nfirst(D2,3)","category":"page"},{"location":"man/Example1/","page":"Example 1","title":"Example 1","text":"function makeprior2(p=-4; μ=[-8.38, -0.22, 0.97*0.27], ip=1.0)\n   Σ = Diagonal([ip*1.37, (10.0^p)*0.02^2, (10.0^p)*0.15^2])\n   return MvNormal(μ, Σ)\nend\n\nfunction lsf(p::Float64, D::DataFrame; lags=0:25)\n    Bprior = makeprior2(p, ip=1000.0)\n    μ = mean(Bout(D, Bprior=Bprior, lags=lags).value[3])\n    return (λ=10.0^p, b₀=μ[1], b₁=μ[2], b₂=μ[3])\nend\n\nDbpar = DataFrame(p=[-4:0.5:2;])\nDb = select(Dbpar, :p=>ByRow(p->lsf(p, D2[:, 2:end], lags=0:100))=>AsTable)\nDb","category":"page"},{"location":"man/Example3/#Example-3","page":"Example 3","title":"Example 3","text":"","category":"section"},{"location":"man/Example3/","page":"Example 3","title":"Example 3","text":"This example shows that the algorithm of Hannart et al. (2014) and the WTLS algorithm in Appendix C and D lead to approximately the same basic ML solution.","category":"page"},{"location":"man/Example3/","page":"Example 3","title":"Example 3","text":"using CSV, DataFrames, OptimizationOptimJL, Distributions\nimport ProfileLikelihood as pl\nimport Measurements as Mm\nusing Measurements: ±\nusing MeasurementErrorModels\n\nMunit =  Mm.Measurement{Float64}\nf1 = joinpath(dirname(pathof(MeasurementErrorModels)), \"../docs/src/assets/\")\nD2 = CSV.read(f1*\"Palaui.csv\", DataFrame, types=Dict([2,3,4].=>Munit))\nfirst(D2, 3)","category":"page"},{"location":"man/Example3/","page":"Example 3","title":"Example 3","text":"\nfunction getci(prof, i)\n    ml = pl.get_likelihood_solution(prof)[i]\n    ci = pl.get_confidence_intervals(prof[i])\n    twoσ = maximum(abs, ml .- [ci.lower, ci.upper])\n    return ml ± 0.5*twoσ\nend\n\nB0 = [-10, -0.22, 0.27]\ndata = dataf(D2[:,2:end], A=\"hannart\")\nprob = pl.LikelihoodProblem(hloglik, B0, data=data, syms=[:b₀, :b₁, :b₂])\nml1 = pl.mle(prob, NelderMead())","category":"page"},{"location":"man/Example3/","page":"Example 3","title":"Example 3","text":"lb, ub = [-50.0, -0.5, 0.0], [0.0, 0.0, 1.0]\nparam_ranges = pl.construct_profile_ranges(ml1, lb, ub, [200, 200, 200])\nprof = pl.profile(prob, ml1; param_ranges=param_ranges, parallel=true)\nBhml = [getci(prof, i) for i in [1 2 3]]\nBml = wtls(D2[:,2:end])\nBwtls = permutedims(mean(Bml) .± sqrt.(var(Bml)))\nDb = DataFrame([Bhml; Bwtls], [:b₀, :b₁, :b₂])\ninsertcols!(Db, 1, :Method=>[\"HML\",\"WTLS\"])","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = MeasurementErrorModels","category":"page"},{"location":"#Measurement-Error-Models","page":"Home","title":"Measurement Error Models","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"MeasurementErrorModels.jl, short for Measurement Error Proxy System Models (MEPSMs), is aimed at calibrating MEPSMs (Fischer 2023). The code could be used generally for multiple-input measurement error models.","category":"page"},{"location":"#[Theory](@ref)","page":"Home","title":"Theory","text":"","category":"section"},{"location":"#[Guide](@ref-Example1)","page":"Home","title":"Guide","text":"","category":"section"},{"location":"#[Library](@ref)","page":"Home","title":"Library","text":"","category":"section"},{"location":"#Reference","page":"Home","title":"Reference","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Fischer, M (2023) Measurement Error Proxy System Models: MEPSM v0.2. link to preprint","category":"page"},{"location":"#Documentation","page":"Home","title":"Documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The documentation was built using Documenter.jl.","category":"page"},{"location":"","page":"Home","title":"Home","text":"println(\"Documentation built with Julia $(VERSION).\") # hide","category":"page"}]
}
