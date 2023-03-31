


# Example 3

This example shows that the algorithm of Hannart et al. (2014) and the WTLS algorithm in Appendix C and D lead to approximately the same basic ML solution.

```@example 3
using CSV, DataFrames, OptimizationOptimJL, Distributions
import ProfileLikelihood as pl
import Measurements as Mm
using Measurements: ±
using MeasurementErrorModels

Munit =  Mm.Measurement{Float64}
f1 = joinpath(dirname(pathof(MeasurementErrorModels)), "../docs/src/assets/")
D2 = CSV.read(f1*"Palaui.csv", DataFrame, types=Dict([2,3,4].=>Munit))
first(D2, 3)
```


```@example 3

function getci(prof, i)
    ml = pl.get_likelihood_solution(prof)[i]
    ci = pl.get_confidence_intervals(prof[i])
    twoσ = maximum(abs, ml .- [ci.lower, ci.upper])
    return ml ± 0.5*twoσ
end

B0 = [-10, -0.22, 0.27]
data = dataf(D2[:,2:end], A="hannart")
prob = pl.LikelihoodProblem(hloglik, B0, data=data, syms=[:b₀, :b₁, :b₂])
ml1 = pl.mle(prob, NelderMead())
```

```@example 3
lb, ub = [-50.0, -0.5, 0.0], [0.0, 0.0, 1.0]
param_ranges = pl.construct_profile_ranges(ml1, lb, ub, [200, 200, 200])
prof = pl.profile(prob, ml1; param_ranges=param_ranges, parallel=true)
Bhml = [getci(prof, i) for i in [1 2 3]]
Bml = wtls(D2[:,2:end])
Bwtls = permutedims(mean(Bml) .± sqrt.(var(Bml)))
Db = DataFrame([Bhml; Bwtls], [:b₀, :b₁, :b₂])
insertcols!(Db, 1, :Method=>["HML","WTLS"])
```
