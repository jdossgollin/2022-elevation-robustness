using LaTeXStrings
using StatsBase

function plot_rcp_tradeoffs(
    u::Matrix{<:HouseElevation.Outcome},
    s::Vector{<:HouseElevation.LSLSim},
    x::Vector{<:Unitful.Length};
    house_value_usd::T,
    house_floor_area::A,
) where {T<:Real,A<:Unitful.Area}

    # define constant
    models = (
        (rcp=2.6, model="BRICK Slow"),
        (rcp=4.5, model="BRICK Fast"),
        (rcp=6.0, model="K14"),
        (rcp=8.5, model="DP16"),
    )

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
        for (model, color) in zip(models, colors)
            w = weights([(si.rcp == model.rcp) & (si.model == model.model) for si in s])
            cond = vec(mean(var, w; dims=1))
            plot!(
                p,
                Δh_ft,
                cond;
                label="RCP $(model.rcp), $(model.model)",
                color=color,
                linewidth=3,
            )
            idx = argmin(cond)
            scatter!(p, [Δh_ft[idx]], [cond[idx]]; label=false, color=color)
        end
        p

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

    savefig(plots_dir("tradeoffs-by-rcp.pdf"))
    return p
end
