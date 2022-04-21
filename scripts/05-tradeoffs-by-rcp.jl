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
    model_rcp_plot = (
        (rcp=2.6, model="BRICK Slow"),
        (rcp=6.0, model="K14"),
        (rcp=8.5, model="DP16"),
    )

    all_rcp = unique([si.rcp for si in s])
    all_models = unique([si.model for si in s])

    total_cost = map(ui -> (ui.led_usd + ui.upfront_cost_usd) / house_value_usd, u)
    Δh_ft = ustrip.(u"ft", x)
    led = map(ui -> ui.led_usd / house_value_usd, u)

    # we need ticks to plot
    x_ticks = 0:2:12 # in feet
    cost_fn = HouseElevation.get_elevation_cost_function()
    prop_cost = cost_fn.(x_ticks .* 1u"ft", house_floor_area) ./ house_value_usd

    # what we're going to save
    p_archive = []
    varnames = [
        "Expected Lifetime Cost [% House Value]",
        "Lifetime Expected Damages [% House Value]",
    ]
    vars = [total_cost, led]

    for (var, varname) in zip(
        vars,
        varnames,
    )
        p = plot(;
            xlabel=L"Height Increase $\Delta h$ [ft]",
            ylabel=varname,
            linewidth=2,
            xticks=(x_ticks, string.(x_ticks)),
            yformatter=pct_formatter,
            top_margin=12.5Plots.mm,
            left_margin=7.5Plots.mm,
            bottom_margin=7.5Plots.mm,
            legend=ifelse(var == first(vars), :topright, false),
        )

        i = 1 # counter for colors
        for rcp in all_rcp
            for model in all_models
                plot_in_color = (rcp=rcp, model=model) in model_rcp_plot
                w = weights([(si.rcp == rcp) & (si.model == model) for si in s])
                cond = vec(mean(var, w; dims=1))
                plot!(
                    p,
                    Δh_ft,
                    cond;
                    label=ifelse(plot_in_color, "RCP $(rcp), $(model)", ""),
                    color=ifelse(plot_in_color, colors[i], :gray),
                    linewidth=ifelse(plot_in_color, 3.5, 0.5),
                )
                if plot_in_color
                    i += 1
                end
            end
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
    add_panel_letters!(p_archive; fontsize=12, loc=(0.1, 0.95))
    p = plot(p_archive...; layout=(1, 2), link=:x, size=(1100, 550))

    savefig(plots_dir("tradeoffs-by-rcp.pdf"))
    return p
end
