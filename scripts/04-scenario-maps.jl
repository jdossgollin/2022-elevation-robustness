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
                :bins => (LinRange(0, maximum(msl_rise_ft), 60), LinRange(0, 4.5, 60)),
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
                push!(
                    kwargs, :xlabel => "LSLR: $syear to $eyear [ft]", :bottommargin => 7.5mm
                )
                kwargs[:xformatter] = x -> x # update it
            end

            # left column has ylabel, right hgas h₀
            if col == 1 # LHS
                push!(
                    kwargs,
                    :ylabel => "Total Cost [% House Value]",
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

    fake_plot = scatter(
        [0, 0],
        [0, 1];
        zcolor=[0, 3],
        xlims=(1, 1.1),
        xshowaxis=false,
        yshowaxis=false,
        colorbar_title="Density of SOWs",
        grid=false,
        clims=clims,
        c=cmap,
        label=false,
        xticks=[],
        yticks=[],
    )

    @show fake_plot[1][:xaxis][:showaxis]
    @show fake_plot[1][:yaxis][:showaxis]
    l = @layout [grid(N_row, N_col) a{0.1w}]

    p = plot(
        vcat([plots[row, col] for row in 1:N_row for col in 1:N_col], fake_plot)...;
        layout=l,
        size=(1500, 1000),
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
    function proportional_led(ui)
        return ui.led_usd / house_value_usd
    end

    plots = []
    for (fn, fn_name) in zip(
        [proportional_lifetime_cost, proportional_led],
        [
            "Expected Total Cost [% House Value]",
            "Lifetime Expected Damages [% House Value]",
        ],
    )
        prop_cost = fn.(u)

        # develop bins to use
        ystep = 0.1
        Δy = ystep / 2
        msl_plot = collect(0.5:0.1:6.0)

        # calculate the indices of the SLR values that go with each bin
        indices = [findall((y .- Δy) .< msl_rise_ft .< (y .+ Δy)) for y in msl_plot]

        # calcualte expected costs for all bins
        expected_cost = collect(
            hcat(
                [
                    [mean(prop_cost[idx, i]) for idx in indices] for (i, xi) in enumerate(x)
                ]...,
            ),
        )

        x_ft = ustrip.(u"ft", x)
        pi = contourf(
            x_ft,
            msl_plot,
            expected_cost;
            ylabel="LSLR: $syear to $eyear [ft]",
            colorbar_title="$fn_name",
            c=cgrad(:plasma; rev=true),
            linewidth=0,
            levels=30,
        )
        push!(plots, pi)
    end
    xlabel!(last(plots), L"Height Increase $\Delta h$ [ft]")

    p = plot(
        plots...; layout=(2, 1), size=(850, 850), link=:x, leftmargin=5mm, rightmargin=5mm
    )
    savefig(p, plots_dir("scenario-map-height-slr.pdf"))
    return p
end
