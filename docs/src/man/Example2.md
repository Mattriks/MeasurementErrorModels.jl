

# Example 2

In this example, the first posterior (`Bpost1` in the table below) "approaches" the basic ML solution (see [Example 3](@ref)), except the prior has ``λ=1.0``. The final posterior (`Bpost2`) includes the full auto- and cross-correlation in ``Σ_A``, so is different from the basic ML solution.


```@example 2
using CSV, DataFrames, Distributions, LinearAlgebra
import Measurements as Mm
using Measurements: ±
using MeasurementErrorModels

Munit =  Mm.Measurement{Float64}
f1 = joinpath(dirname(pathof(MeasurementErrorModels)), "../docs/src/assets/")
D2 = CSV.read(f1*"Palaui.csv", DataFrame, types=Dict([2,3,4].=>Munit))
first(D2, 3)
```

```@example 2

function marginals(D::AbstractVector)
    z = mean(D[1]) .± sqrt.(var(D[1]))
    return (b₀=z[1], b₁=z[2], b₂=z[3])
end

μ₀, Σ₀ = [-0.22, 0.97*0.27], Diagonal([0.02^2, 0.15^2])
Bprior = Priors.makeprior(D2[:,2:4], MvNormal(μ₀, Σ₀), ip=1000.0)

Dbs = Bout(D2[:,2:end], Bprior=Bprior, lags=0:100)
combine(groupby(Dbs, :Dist), :value=>marginals=>AsTable)
```

