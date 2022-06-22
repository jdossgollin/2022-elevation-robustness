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
            p, 0:0.05:40, prior.dist; label="$(prior.rt) Year Surge", linewidth=3, color=c
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
        ylabel="Ann-Max Storm Surge [ft]",
        left_margin=5mm,
        markerstrokewidth=0,
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
        color = storm.is_tc ? colors[1] : colors[2]
        plot!([x0, x1], [y0, y1]; color=color, linewidth=1.5, label="", markerstrokewidth=0)
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
        prior_fits, plots_dir("surge-prior-mcmc-diagnostics.tex")
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

    # plot the quantiles in gray
    p = plot_return_period(post_gevs)

    # add true line
    rts = range(1.25, 275; length=250) # return periods
    aeps = 1 .- 1 ./ rts # annual exceedance probability
    plot!(p, rts, quantile.(dist, aeps); linewidth=2, color=:blue, label="True")

    # add data
    xp, ys = HouseElevation.weibull_plot_pos(fake_data)
    scatter!(p, 1 ./ xp, ys; label="Observations", color=:gray, alpha=1)

    plot!(p; legend=:topleft)

    # save
    savefig(p, plots_dir("surge-synthetic-data-experiment.pdf"))
    return p
end

# Plot the posterior draws!
function plot_surge_posterior_chains(fits::T) where {T<:HouseElevation.MCMCChains.Chains}
    display(fits)

    HouseElevation.write_diagnostics(
        fits, plots_dir("surge-posterior-mcmc-diagnostics.tex")
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

    # add gray lines
    p = plot_return_period(post_gevs)
    y_ticks = [5, 7.5, 10, 12.5, 15]
    plot!(p; yaxis=:log, yticks=(y_ticks, string.(y_ticks)))

    # add weibul plot pls
    xp, ys = HouseElevation.weibull_plot_pos(surge_ft)
    scatter!(p, 1 ./ xp, ys; label="Observations", color=:gray, alpha=1)
    return p
end

# make a plot of the return period
function plot_surge_prior_return() where {T<:HouseElevation.MCMCChains.Chains}
    prior_model = HouseElevation.StationaryGEV(missing)
    prior_fits = HouseElevation.get_posterior(
        prior_model, "prior_model", 100_000; n_chains=4, overwrite=false, drop_warmup=true
    )

    # fit the model
    df = DataFrames.DataFrame(prior_fits)
    gevs = [GeneralizedExtremeValue(row[:μ], row[:σ], row[:ξ]) for row in eachrow(df)]

    p = plot_return_period(gevs; type="Prior")
    savefig(p, plots_dir("surge-prior-return.pdf"))
    return p
end

function plot_surge_obs_return(
    annual::HouseElevation.AnnualGageRecord, fits::T
) where {T<:HouseElevation.MCMCChains.Chains}
    p1 = plot_annmax_floods(annual)
    xlabel!(p1, "Time [year]")
    p2 = plot_surge_posterior_return(annual, fits)
    p2 = plot(p2; ylabel="")

    y_ticks = collect(3:10)
    p1 = plot!(p1; yaxis=:log, yticks=(y_ticks, string.(y_ticks)))
    p2 = plot!(p2; yaxis=:log, yticks=(y_ticks, string.(y_ticks)))

    fontsize = 10
    annotate!(p1, (0.05, 0.95), text("(A)", :left, fontsize))
    annotate!(p2, (1/1000, 0.95), text("(B)", :left, fontsize))
    p = plot(
        p1,
        p2;
        link=:y,
        ylims=(2.5, 10),
        layout=grid(1, 2; widths=[0.6, 0.35]),
        size=[1000, 350],
        bottommargin=8mm,
        leftmargin=8mm,
    )
    savefig(p, plots_dir("surge-obs-return.pdf"))
    return p
end
