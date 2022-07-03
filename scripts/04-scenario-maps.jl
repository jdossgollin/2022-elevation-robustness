using ColorSchemes
using LaTeXStrings
using Plots
using Plots: mm
using Unitful

function plot_scenario_map_slr_cost(;
    x_plot::Vector{<:L},
    elevation_init_bfe::Vector{<:L},
    fits::C,
    house_value_usd::T,
    syear::Int,
    eyear::Int,
    house_floor_area::A,
    discount_rate::T,
    overwrite::Bool=false,
) where {L<:Unitful.Length,A<:Unitful.Area,T<:Real,C<:HouseElevation.MCMCChains.Chains}

    # get the SLR
    s = HouseElevation.get_lsl(; syear=syear, eyear=eyear)

    # calc the base flood elevation
    bfe = calc_bfe(fits, s, syear) # has 1% chance of flood in 2022

    # on the x axis we will plot MSL from `syear`` to `eyear``
    msl_rise = get_year_data(s, eyear)
    msl_rise_ft = ustrip.(u"ft", msl_rise)

    N_row = length(elevation_init_bfe)
    N_col = length(x_plot)
    plots = hcat([[plot() for row in 1:N_row] for col in 1:N_col]...)

    clims = (0, 8_000)
    cmap = cgrad(:plasma; scale=exp, rev=true)

    for row in 1:N_row
        h0_bfe = elevation_init_bfe[row]
        h0 = h0_bfe + bfe
        for col in 1:N_col
            Δh = x_plot[col]
            p = plots[row, col]
            u = HouseElevation.get_outcomes(
                [Δh];
                elevation_init=h0,
                syear=syear,
                eyear=eyear,
                fits=fits,
                house_floor_area=house_floor_area,
                discount_rate=discount_rate,
                house_value_usd=house_value_usd,
                overwrite=overwrite,
            )
            total_cost_usd = [ui.led_usd + ui.upfront_cost_usd for ui in u][:, 1]
            total_cost_prop = total_cost_usd ./ house_value_usd

            # a way to add stuff based on where in the plot we are
            kwargs = Dict{Symbol,Any}(
                :label => false,
                :xformatter => x -> "",
                :bins => (LinRange(0, maximum(msl_rise_ft), 40), LinRange(0, 4.5, 40)),
                :c => cmap,
                :clims => clims,
                :colorbar => false,
            )

            # top row has the titles and the bottom has xlabels
            if row == 1
                push!(
                    kwargs,
                    :title => L"$\Delta h$ = %$Δh",
                    :topmargin => 5mm,
                    :title_align => :center,
                )
            elseif row == N_row
                kwargs[:xformatter] = identity
                push!(
                    kwargs,  :bottommargin => 7.5mm
                )
                if col == 2
                    push!(
                    kwargs, :xlabel => "SLR: $syear to $eyear [ft]"
                )
                end
            end

            # left column has ylabel, right hgas h₀
            if col == 1 # LHS
                push!(
                    kwargs,
                    :ylabel => "Total Cost\n[% House Value]",
                    :leftmargin => 7.5mm,
                    :yformatter => pct_formatter,
                )
            elseif col == N_col # RHS
                ft = u"ft"
                if h0_bfe > 0u"ft"
                    ylabel = L"$h_0$ = %$(round(ustrip(ft, h0_bfe), digits=1)) ft above BFE"
                else
                    ylabel = L"$h_0$ = %$(round(ustrip(ft, -h0_bfe), digits=1)) ft below BFE"
                end
                push!(
                    kwargs,
                    :ylabel => ylabel,
                    :guide_position => :right,
                    :yguidefontrotation => 180,
                    :rightmargin => 12.5mm,
                    :yformatter => y -> "",
                )
            else
                push!(kwargs, :yformatter => y -> "")
            end

            #scatter!(p, msl_rise_ft, total_cost_prop; kwargs...)
            histogram2d!(p, msl_rise_ft, total_cost_prop; kwargs...)
        end
    end

    cbar_fake_plot = scatter(
        [0, 0],
        [0, 1];
        zcolor=[0, 3],
        xlims=(1, 1.1),
        xshowaxis=false,
        yshowaxis=false,
        colorbar_title="Count of SOWs",
        colorbar_titlefontrotation=90,
        grid=false,
        clims=clims,
        c=cmap,
        label=false,
        xticks=[],
        yticks=[],
    )

    cbar_fake_plot[1][:xaxis][:showaxis]
    cbar_fake_plot[1][:yaxis][:showaxis]

    harrow = plot(
        [0.2, 0.9],
        [0, 0];
        axis=false,
        ticks=false,
        legend=false,
        arrow=:closed,
        color=:black,
        xlims=(0, 1),
        ylims=(-0.05, 0.25),
    )
    annotate!(
        harrow, [0.5], [0.1], text(L"Increasing House Elevation ($\Delta h$)", :center, 14)
    )

    varrow = plot(
        [0, 0],
        [0.1, 0.9];
        axis=false,
        ticks=false,
        legend=false,
        arrow=:closed,
        color=:black,
        xlims=(-0.25, 0.05),
        ylims=(0, 1),
    )
    annotate!(
        varrow,
        [-0.1],
        [0.5],
        text(L"Increasing Initial Exposure ($h_0$)", :center, 14; rotation=90),
    )

    l = @layout [
        a{0.075h}
        [b{0.05w} grid(N_row, N_col) c{0.125w}]
    ]

    add_panel_letters!(plots)
    p = plot(
        vcat(
            harrow,
            varrow,
            [plots[row, col] for row in 1:N_row for col in 1:N_col],
            cbar_fake_plot,
        )...;
        layout=l,
        size=(1200, 750),
        dpi=250,
        link=:y,
    )
    savefig(p, plots_dir("scenario-map-slr-cost.pdf"))
    return p
end

function plot_scenario_map_height_slr(;
    u::Array{<:HouseElevation.Outcome},
    s::Vector{<:HouseElevation.LSLSim},
    x::Vector{<:Unitful.Length},
    house_value_usd,
)
    # calculate first and last years
    s1 = first(s)
    syear = minimum(s1.years)
    eyear = maximum(s1.years)

    # get the sea level data
    msl_rise = get_year_data(s, eyear)
    msl_rise_ft = ustrip.(u"ft", msl_rise)

    function proportional_lifetime_cost(ui)
        return (ui.led_usd + ui.upfront_cost_usd) / house_value_usd
    end

    x_ticks = 0:2:10 # in feet

    prop_cost = proportional_lifetime_cost.(u)

    # develop bins to use
    ystep = 0.1
    Δy = ystep / 2
    msl_plot = collect(0.5:0.1:9.5)

    # calculate the indices of the SLR values that go with each bin
    indices = [findall((y .- Δy) .< msl_rise_ft .< (y .+ Δy)) for y in msl_plot]

    # calcualte expected costs for all bins
    expected_cost = collect(
        hcat([[mean(prop_cost[idx, i]) for idx in indices] for (i, xi) in enumerate(x)]...)
    )

    x_ft = ustrip.(u"ft", x)

    p = heatmap(
        x_ft,
        msl_plot,
        expected_cost;
        ylabel="SLR from $syear to $eyear [ft]",
        colorbar_title="Expected Total Costs\n[% House Value]",
        colorbar_formatter=pct_formatter,
        c=cgrad(:plasma; rev=true),
        linewidth=0,
        levels=LinRange(0.25, 3.25, 31),
        #clim=(0.25, 3.25),
        xticks=(x_ticks, string.(x_ticks)),
        xlabel=L"Height Increase $\Delta h$ [ft]",
        size=(800, 450),
        bottom_margin=5mm,
        left_margin=5mm,
        right_margin=7.5mm,
    )
    contour!(p, x_ft, msl_plot, expected_cost; linewidth=0.5, levels=25)
    savefig(p, plots_dir("scenario-map-height-slr.pdf"))
    return p
end
