using CSV
using DataFrames
using Distributions
using LaTeXStrings
using Plots
using Plots: mm
using Random
using StatsBase
using StatsPlots
using Unitful

# make a plot of all the GEV priors and save this plot
function plot_surge_gev_priors()
    # first, set up the plot
    p = plot(;
        xlabel="Storm Surge [ft]",
        ylabel="Prior Predictive Density",
        leftmargin=5mm,
        legend=:topright,
        yticks=0:0.2:1,
    )
    # loop through all the distributions
    for (prior, c) in zip(HouseElevation.gev_priors, colors)
        plot!(
            p, 0:0.025:25, prior.dist; label="$(prior.rt) Year Surge", linewidth=3, color=c
        )
    end
    savefig(p, plots_dir("surge-gev-priors.pdf"))
    return p
end

# plot the historic annual maximum floods
function plot_annmax_floods(annual::HouseElevation.AnnualGageRecord)

    # define some historic storms
    historic_norfolk_storms = [
        (name="Chesapeake-Potomac", year=1933, is_tc=true, Δx=7.5, Δy=0.75),
        (name="Outer Banks hurricane", year=1936, is_tc=true, Δx=10, Δy=0.75),
        (name="Isabel", year=2003, is_tc=true, Δx=-7.5, Δy=0.5),
        (name="Irene", year=2015, is_tc=true, Δx=2.5, Δy=1),
        (name="Nor'Ida", year=2009, is_tc=false, Δx=0, Δy=1),
        (name="Ash Wednesday", year=1962, is_tc=false, Δx=10, Δy=0.75),
    ]

    # Sewells Point, VA is default

    surge_ft = ustrip.(u"ft", annual.max_surge) # scalarize in ft

    # plot the historical storms
    p = scatter(
        annual.t,
        surge_ft;
        label=false,
        markercolor=:gray,
        ylabel="Ann-Max Storm Surge at $(annual.stn.gage_name) [ft]",
        left_margin=5mm,
    )

    # show historic storms
    for storm in historic_norfolk_storms
        yobs = ustrip(surge_ft[findfirst(annual.t .== storm.year)])
        x_text = storm.year + storm.Δx
        y_text = yobs + storm.Δy
        x0 = x_text + 0.2 * (storm.year - x_text)
        x1 = x_text + 0.9 * (storm.year - x_text)
        y0 = y_text + 0.2 * (yobs - y_text)
        y1 = y_text + 0.9 * (yobs - y_text)
        color = storm.is_tc ? colors[5] : colors[6]
        plot!([x0, x1], [y0, y1]; color=color, linewidth=1.5, label="")
        annotate!(p, x_text, y_text, text(storm.name, :center, 7; color=color))
    end
    return p
end

# plot the samples from the prior
function plot_surge_prior_chains()
    prior_model = HouseElevation.StationaryGEV(missing)
    prior_fits = HouseElevation.get_posterior(
        prior_model, "prior_model", 100_000; n_chains=4, overwrite=false, drop_warmup=true
    )
    display(prior_fits)

    HouseElevation.write_diagnostics(
        prior_fits, plots_dir("surge-prior-mcmc-diagnostics.csv")
    )

    p = plot(prior_fits)
    savefig(p, plots_dir("surge-prior-chains.png"))
    return p
end

# Plot the 
function plot_surge_synthetic_experiment(annual::HouseElevation.AnnualGageRecord)

    # create fake data
    function make_fake_data(dist, N)
        rng = Random.MersenneTwister(713)
        return rand(rng, dist, N)
    end
    N = length(annual.t)

    # get fake data
    dist = Distributions.GeneralizedExtremeValue(4, 0.5, 0.15)
    fake_data = make_fake_data(dist, N)

    # fit the model
    posterior = HouseElevation.get_posterior(
        HouseElevation.StationaryGEV(fake_data), "fake_data", 25_000
    )
    posterior_df = DataFrames.DataFrame(posterior)
    post_gevs = [
        GeneralizedExtremeValue(row[:μ], row[:σ], row[:ξ]) for row in eachrow(posterior_df)
    ]

    # set up the return period plot
    rts = range(1.25, 300; length=250) # return periods
    aeps = 1 .- 1 ./ rts # annual exceedance probability
    xticks = [2, 5, 10, 25, 50, 100, 250]

    ub1 = [quantile([quantile(d, xi) for d in post_gevs], 0.95) for xi in aeps]
    lb1 = [quantile([quantile(d, xi) for d in post_gevs], 0.05) for xi in aeps]
    ub2 = [quantile([quantile(d, xi) for d in post_gevs], 0.9) for xi in aeps]
    lb2 = [quantile([quantile(d, xi) for d in post_gevs], 0.1) for xi in aeps]
    ub3 = [quantile([quantile(d, xi) for d in post_gevs], 0.75) for xi in aeps]
    lb3 = [quantile([quantile(d, xi) for d in post_gevs], 0.25) for xi in aeps]

    p = plot(;
        xlabel="Return Period [years]",
        ylabel="Return Level [ft]",
        xscale=:log,
        legend=:topleft,
        xticks=(xticks, string.(xticks)),
        dpi=250, # for saving
        size=(500, 500),
        title="Fake Data Experiment: N=$N",
    )
    plot!(
        rts,
        ub1;
        fillbetween=lb1,
        fillcolor=:gray,
        fillalpha=0.35,
        linecolor=false,
        label="95% Posterior CI",
    )
    plot!(
        rts,
        ub2;
        fillbetween=lb2,
        fillcolor=:gray,
        fillalpha=0.35,
        linecolor=false,
        label="80% Posterior CI",
    )
    plot!(
        rts,
        ub3;
        fillbetween=lb3,
        fillcolor=:gray,
        fillalpha=0.35,
        linecolor=false,
        label="50% Posterior CI",
    )
    plot!(p, rts, quantile.(dist, aeps); linewidth=2, color=:blue, label="True")
    savefig(p, plots_dir("surge-synthetic-data-experiment.pdf"))
    return p
end

# Plot the posterior draws!
function plot_surge_posterior_chains(fits::T) where {T<:HouseElevation.MCMCChains.Chains}
    display(fits)

    HouseElevation.write_diagnostics(
        fits, plots_dir("surge-posterior-mcmc-diagnostics.csv")
    )

    p = plot(fits)
    savefig(p, plots_dir("surge-posterior-chains.png"))
    return p
end

function plot_surge_posterior_teststats(
    annual::HouseElevation.AnnualGageRecord, fits::T
) where {T<:HouseElevation.MCMCChains.Chains}

    # function to plot a test statistics
    function plot_test_stat(t, y, yhat; title="", xlabel="", draw_label::Bool=false)
        t_posterior = [t(yi) for yi in yhat]
        observed = t(y)
        bins = range(quantile(t_posterior, 0.002), quantile(t_posterior, 0.998); length=50)
        bin_label = draw_label ? L"PPD: $p(\tilde{y} | y)$" : false
        line_label = draw_label ? L"Obs: $y$" : false
        p = plot(;
            title_align=:left,
            yticks=:none,
            yaxis=([], false),
            title=title,
            xlabel=xlabel,
            titlefontsize=12,
        )
        histogram!(
            p,
            t_posterior;
            bins=bins,
            label=bin_label,
            normalize=:pdf,
            fillcolor=false,
            linecolor=:black,
            linewidth=0.5,
        )
        vline!(p, [observed]; linewidth=4, label=line_label, color=colors[1])
        return p
    end
    function MannKendall(x::Vector{<:Real})
        return StatsBase.corkendall([1:length(x);], x)
    end
    function sample_predictive(fits, N)
        return rand.(
            GeneralizedExtremeValue.(vec(fits[:μ]), vec(fits[:σ]), vec(fits[:ξ])), N
        )
    end

    # get the data and the fake data
    surge_ft = ustrip.(u"ft", annual.max_surge) # scalarize in ft
    N = length(annual.t)
    ŷ = sample_predictive(fits, N)

    # what test stats will we plot?
    test_stats = [
        ("(a) Lag 1 PACF", "Autocorrelation", x -> StatsBase.pacf(x, 1:1)[1]),
        ("(b) Lag 2 PACF", "Autocorrelation", x -> StatsBase.pacf(x, 2:2)[1]),
        ("(c) Largest Flood (in $N Years)", "Surge [ft]", maximum),
        ("(d) Smallest Flood (in $N Years)", "Surge [ft]", minimum),
        ("(e) Median Flood", "Surge [ft]", median),
        ("(f) Mann-Kendall Test", "Test Statistic", MannKendall),
    ]

    # make the plots
    test_stat_plots = []
    for (i, test_stat) in enumerate(test_stats)
        title, xlabel, t = test_stat
        pi = plot_test_stat(t, surge_ft, ŷ; title=title, xlabel=xlabel, draw_label=i == 1)
        push!(test_stat_plots, pi)
    end
    p = plot(test_stat_plots...; size=(900, 600))
    savefig(p, plots_dir("surge-test-statistics.pdf"))
    return p
end

# make a plot of the return period
function plot_surge_posterior_return(
    annual::HouseElevation.AnnualGageRecord, fits::T
) where {T<:HouseElevation.MCMCChains.Chains}

    # get the raw data
    surge_ft = ustrip.(u"ft", annual.max_surge)
    N = length(surge_ft)

    # fit the model
    posterior_df = DataFrames.DataFrame(fits)
    post_gevs = [
        GeneralizedExtremeValue(row[:μ], row[:σ], row[:ξ]) for row in eachrow(posterior_df)
    ]

    rts = range(1.25, 275; length=250) # return periods
    aeps = 1 .- 1 ./ rts # annual exceedance probability
    xticks = [2, 5, 10, 25, 50, 100, 250]

    ub1 = [quantile([quantile(d, xi) for d in post_gevs], 0.95) for xi in aeps]
    lb1 = [quantile([quantile(d, xi) for d in post_gevs], 0.05) for xi in aeps]
    ub2 = [quantile([quantile(d, xi) for d in post_gevs], 0.9) for xi in aeps]
    lb2 = [quantile([quantile(d, xi) for d in post_gevs], 0.1) for xi in aeps]
    ub3 = [quantile([quantile(d, xi) for d in post_gevs], 0.75) for xi in aeps]
    lb3 = [quantile([quantile(d, xi) for d in post_gevs], 0.25) for xi in aeps]

    p = plot(;
        xlabel="Return Period [years]",
        ylabel="Return Level [ft]",
        xscale=:log,
        legend=:bottomright,
        xticks=(xticks, string.(xticks)),
    )
    plot!(
        rts,
        ub1;
        fillbetween=lb1,
        fillcolor=:gray,
        fillalpha=0.35,
        linecolor=false,
        label="95% Posterior CI",
    )
    plot!(
        rts,
        ub2;
        fillbetween=lb2,
        fillcolor=:gray,
        fillalpha=0.35,
        linecolor=false,
        label="80% Posterior CI",
    )
    plot!(
        rts,
        ub3;
        fillbetween=lb3,
        fillcolor=:gray,
        fillalpha=0.35,
        linecolor=false,
        label="50% Posterior CI",
    )

    xp, ys = HouseElevation.weibull_plot_pos(surge_ft)
    scatter!(p, 1 ./ xp, ys; label="Obs (Weibull Plot Pos.)", color=colors[1], alpha=1)
    savefig(p, plots_dir("surge-synthetic-data-experiment.pdf"))
    return p
end

function plot_surge_obs_return(
    annual::HouseElevation.AnnualGageRecord, fits::T
) where {T<:HouseElevation.MCMCChains.Chains}
    p1 = plot_annmax_floods(annual)
    xlabel!(p1, "Time [year]")
    p2 = plot_surge_posterior_return(annual, fits)
    p2 = plot(p2; ylabel="")
    add_panel_letters!([p1]; fontsize=14)
    annotate!(p2, [1], [8.5], text("(B)", :left, 14))
    p = plot(
        p1,
        p2;
        link=:y,
        ylims=(2, 9),
        layout=grid(1, 2; widths=[0.6, 0.35]),
        size=[1000, 350] .* 1.35,
        bottommargin=8mm,
        leftmargin=8mm,
    )
    savefig(p, plots_dir("surge-obs-return.pdf"))
    return p
end

function make_surge_plots(
    annual::HouseElevation.AnnualGageRecord, fits::T
) where {T<:HouseElevation.MCMCChains.Chains}
    plot_surge_gev_priors()
    plot_annmax_floods(annual)
    plot_surge_prior_chains()
    plot_surge_synthetic_experiment(annual)
    plot_surge_posterior_chains(fits)
    plot_surge_posterior_teststats(annual, fits)
    plot_surge_posterior_return(annual, fits)
    plot_surge_obs_return(annual, fits)
    return nothing
end
