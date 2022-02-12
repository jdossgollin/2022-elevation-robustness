using LaTeXStrings
using Plots

function plot_grid_scheme()
    dist = Distributions.Normal()
    ψ = quantile.(dist, [0.05, 0.3, 0.7, 0.95])
    ψ_labels = [L"$\psi_%$i$" for i in 1:length(ψ)]
    xlimits = (-3, 3)

    ymax = pdf(dist, 0)
    grids = [0.5 * (ψ[j] + ψ[j + 1]) for j in 1:(length(ψ) - 1)]
    grid_labels = [L"$\frac{\psi_%$j + \psi_%$(j+1)}{2}$" for j in 1:length(grids)]

    p = plot(
        x -> cdf(dist, x);
        xticks=(vcat(ψ, grids), vcat(ψ_labels, grid_labels)),
        linewidth=3,
        color=colors[1],
        label="True Distribution",
        xlims=xlimits,
        ylabel="CDF",
        xlabel="Mean Sea Level in 2100 [ft]",
        legend=:right,
        bottom_margin=5Plots.mm,
        left_margin=5Plots.mm,
        size=(800, 350),
    )
    scatter!(p, ψ, zeros(size(ψ)); label="Observations", color=colors[2], markersize=5)

    ymax_prev = 0.0
    x_arrow = range(minimum(xlimits), 0; length=length(ψ) + 2)[2:(end - 1)]
    for i in 1:length(ψ)
        xmin = minimum(xlimits)
        if i == length(ψ)
            xmax = maximum(xlimits)
        else
            xmax = grids[i]
        end
        ymax = cdf(dist, xmax)
        plot!(
            p,
            [minimum(xlimits), xmax],
            [ymax, ymax];
            label=false,
            color=colors[6],
            linestyle=:dot,
        )
        plot!(p, [xmax, xmax], [0, ymax]; label=false, color=colors[3], linestyle=:dash)
        plot!(
            p,
            [x_arrow[i], x_arrow[i]],
            [ymax_prev, ymax];
            color=colors[5],
            linewidth=3,
            label=false,
        )
        annotate!(
            p,
            [x_arrow[i] + 0.2],
            [(ymax_prev + ymax) / 2],
            text(L"$w_%$i$", :center, 10; color=colors[5]),
        )
        ymax_prev = ymax
    end
    savefig(p, plots_dir("grid-sketch.pdf"))
    return p
end

function get_priors()
    priors = [
        (name="Slow SLR", dist=Gamma(2, 0.5)),
        (name="Intermediate", dist=Gamma(1.5, 1.5)),
        (name="Fast SLR", dist=Gamma(4.5, 1.25)),
    ]
    return priors
end

function plot_priors()
    priors = get_priors()
    p = plot(;
        xlabel="Sea Level Rise, 2022-2100, at Sewells Point, VA [ft]",
        ylabel="Probability Density",
    )
    for (prior, color) in zip(priors, colors)
        plot!(p, prior.dist, 0, 12.5; label=prior.name, color=color, linewidth=3)
    end
    p
    savefig(p, plots_dir("lsl-priors.pdf"))
    return p
end

function plot_prior_tradeoffs(
    u::Matrix{<:HouseElevation.Outcome},
    s::Vector{<:HouseElevation.LSLSim},
    x::Vector{<:Unitful.Length};
    house_value_usd::T,
    house_floor_area::A,
) where {T<:Real,A<:Unitful.Area}

    # define priors over SLR in 2100, in ft
    priors = get_priors()

    total_cost = map(ui -> (ui.led_usd + ui.upfront_cost_usd) / house_value_usd, u)
    Δh_ft = ustrip.(u"ft", x)
    led = map(ui -> ui.led_usd / house_value_usd, u)

    # we need ticks to plot
    x_ticks = 0:2:14 # in feet
    cost_fn = HouseElevation.get_elevation_cost_function()
    prop_cost = cost_fn.(x_ticks .* 1u"ft", house_floor_area) ./ house_value_usd

    # what we're going to save
    p_archive = []

    for (var, varname) in zip(
        [total_cost, led],
        [
            "Expected Lifetime Cost [% House Value]",
            "Lifetime Expected Damages [% House Value]",
        ],
    )
        p = plot(;
            xlabel=L"Height Increase $\Delta h$ [ft]",
            ylabel=varname,
            linewidth=2,
            xticks=(x_ticks, string.(x_ticks)),
            yformatter=y -> pct_formatter(y),
            top_margin=12.5Plots.mm,
            left_margin=7.5Plots.mm,
            bottom_margin=7.5Plots.mm,
            legend=:bottomright,
        )
        for (prior, color) in zip(priors, colors)
            w = HouseElevation.make_weights(s, prior.dist)
            cond = vec(mean(var, w; dims=1))
            plot!(p, Δh_ft, cond; label=prior.name, color=color, linewidth=3)
            idx = argmin(cond)
            scatter!(p, [Δh_ft[idx]], [cond[idx]]; label=false, color=color)
        end

        # add the cost on upper x axis
        p = plot!(
            twiny(p),
            Δh_ft,
            [mean(var[:, i]) for (i, _) in enumerate(Δh_ft)];
            xticks=(x_ticks, string.(round.(prop_cost, digits=2))),
            xlabel="Up-Front Cost [% House Value]",
            linewidth=0,
            alpha=0,
            label=false,
            yticks=false,
        )

        # add to our plots
        push!(p_archive, p)
    end
    add_panel_letters!(p_archive; fontsize=12)
    p = plot(p_archive...; layout=(1, 2), link=:x, size=(1200, 600))

    savefig(plots_dir("tradeoffs-by-prior.pdf"))
    return p
end

"""
Plot the weights
"""
function plot_weight(s::Vector{<:HouseElevation.LSLSim})
    # define priors over SLR in 2100, in ft
    priors = get_priors()
    p = plot()
    df = DataFrame(:model => ["RCP $(si.rcp), $(si.model)" for si in s])
    for prior in priors
        w = HouseElevation.make_weights(s, prior.dist)
        df[!, prior.name] = w
    end
    df2 = stack(df, Not([:model]); variable_name=:prior, value_name=:weight)
    w_avg = combine(groupby(df2, [:model, :prior]), :weight => mean => :weight)
    p = @df w_avg groupedbar(
        :model,
        :weight,
        group=:prior,
        palette=colors,
        legend=:topright,
        legendtitle="Prior",
        yformatter=blank_formatter,
        ylabel="Average Weight Applied",
        xrotation=30,
        bottom_margin=10Plots.mm,
        left_margin=10Plots.mm,
    )
    return p
end
