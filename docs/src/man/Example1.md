
# [Example 1](@id Example1)

In this example, the final posterior is calculated for a range of ``λ`` values.  
The solution for ``λ=1.0`` is the same as `Bpost2` in [Example 2](@ref).

```@example 1

using CSV, DataFrames, Distributions, LinearAlgebra
import Measurements as Mm
using MeasurementErrorModels

Munit =  Mm.Measurement{Float64}
f1 = joinpath(dirname(pathof(MeasurementErrorModels)), "../docs/src/assets/")
D2 = CSV.read(f1*"Palaui.csv", DataFrame, types=Dict([2,3,4].=>Munit))
first(D2,3)
```


```@example 1
function makeprior2(p=-4; μ=[-8.38, -0.22, 0.97*0.27], ip=1.0)
   Σ = Diagonal([ip*1.37, (10.0^p)*0.02^2, (10.0^p)*0.15^2])
   return MvNormal(μ, Σ)
end

function lsf(p::Float64, D::DataFrame; lags=0:25)
    Bprior = makeprior2(p, ip=1000.0)
    μ = mean(Bout(D, Bprior=Bprior, lags=lags).value[3])
    return (λ=10.0^p, b₀=μ[1], b₁=μ[2], b₂=μ[3])
end

Dbpar = DataFrame(p=[-4:0.5:2;])
Db = select(Dbpar, :p=>ByRow(p->lsf(p, D2[:, 2:end], lags=0:100))=>AsTable)
```

