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
    s = HouseElevation.get_norfolk_brick(; syear=syear, eyear=eyear)

    # get the amount of SLR over the time domain
    s1 = first(s)
    syear = minimum(s1.years)
    eyear = maximum(s1.years)

    # calc the base flood elevation
    bfe = calc_bfe(fits, s, syear) # has 1% chance of flood in 2022

    # on the x axis we will plot MSL from `syear`` to `eyear``
    msl_rise = get_year_data(s, eyear) .- get_year_data(s, syear)
    msl_rise_ft = ustrip.(u"ft", msl_rise)

    N_row = length(elevation_init_bfe)
    N_col = length(x_plot)
    plots = hcat([[plot() for row in 1:N_row] for col in 1:N_col]...)

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
            total_cost_prop = total_cost_usd ./ house_value_usd .* 100

            # a way to add stuff based on where in the plot we are
            kwargs = Dict{Symbol,Any}(
                :label => false,
                :markercolor => :black,
                :markersize => 0.5,
                :xformatter => x -> "",
            )

            # top row has the titles and the bottom has xlabels
            if row == 1
                push!(
                    kwargs, :title => "Δh = $Δh", :topmargin => 5mm, :title_align => :center
                )
            elseif row == N_row
                push!(
                    kwargs, :xlabel => "LSLR: $syear to $eyear [ft]", :bottommargin => 7.5mm
                )
                kwargs[:xformatter] = x -> x # update it
            end

            # left column has ylabel, right hgas h₀
            if col == 1 # LHS
                push!(kwargs, :ylabel => "Total Cost [% House Value]", :leftmargin => 7.5mm)
                # custom ticks
                # y_ticks = [10, 25, 50, 100, 250, 500, 1000]
                # kwargs[:yticks] = (y_ticks, string.(y_ticks))
            elseif col == N_col # RHS
                if h0_bfe > 0u"ft"
                    ylabel = "h₀ = $(round(ustrip(u"ft", h0_bfe), digits=1)) ft above BFE"
                else
                    ylabel = "h₀ = $(round(ustrip(u"ft", -h0_bfe), digits=1)) ft below BFE"
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

            scatter!(p, msl_rise_ft, total_cost_prop; kwargs...)
        end
    end
    p = plot(
        [plots[row, col] for row in 1:N_row for col in 1:N_col]...;
        layout=(N_row, N_col),
        size=(1500, 1000),
        link=:y,
    )
    savefig(p, plots_dir("scenario-map-slr-cost.png"))
    return p
end
