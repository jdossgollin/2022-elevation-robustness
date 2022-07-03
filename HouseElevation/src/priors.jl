using DataFrames
using Distributions
using Latexify

"""
Get the ubjective priors for sea level rise

The values give SLR from 2022 to 2100, in ft
"""
function get_priors()
    priors = [
        (name="Slow SLR", dist=Gamma(1.75, 0.5)),
        (name="Uncertain SLR", dist=Gamma(1.75, 1.25)),
        (name="Rapid SLR", dist=Gamma(3.5, 1.25)),
    ]
    return priors
end

"""Write the priors to file"""
function write_priors(fname::String)
    priors = get_priors()
    df = vcat(
        [
            DataFrame([
                "Name" => String(p.name),
                L"$\alpha$" => p.dist.α,
                L"$\theta$" => p.dist.θ,
                [
                    "Q$(q*100)" => round(quantile(p.dist, q); digits=2) for
                    q in [0.025, 0.25, 0.5, 0.75, 0.975]
                ]...,
            ],) for p in priors
        ]...,
    )
    tex_str = Latexify.latexify(
        df; env=:table, fmt="%.2f", booktabs=true, adjustment=:l, latex=false
    )
    open(fname, "w") do io
        write(io, tex_str)
    end
end
