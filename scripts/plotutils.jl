using Colors
using Distributions
using Plots

function add_panel_letters!(panels; fontsize=9, loc=(0.05, 0.95))
    letters = 'A':'Z'
    for (p, letter) in zip(panels, letters)
        annotate!(p, loc, text("($letter)", :left, fontsize))
    end
end

function calc_bfe(fits, all_trajs, start_year)
    dist = Distributions.MixtureModel([
        Distributions.GeneralizedExtremeValue(μ, σ, ξ) for
        (μ, σ, ξ) in zip(vec(fits[:μ]), vec(fits[:σ]), vec(fits[:ξ]))
    ])
    surge_100_year = quantile(dist, 0.99)u"ft"
    lsl_present = mean(get_year_data(all_trajs, start_year))
    return surge_100_year + lsl_present
end

# https://stackoverflow.com/questions/64176617/julia-two-x-axes-for-plot-of-same-data
function twiny(sp::Plots.Subplot)
    sp[:top_margin] = max(sp[:top_margin], 30Plots.px)
    plot!(sp.plt; inset=(sp[:subplot_index], bbox(0, 0, 1, 1)))
    twinsp = sp.plt.subplots[end]
    twinsp[:xaxis][:mirror] = true
    twinsp[:background_color_inside] = RGBA{Float64}(0, 0, 0, 0)
    Plots.link_axes!(sp[:yaxis], twinsp[:yaxis])
    return twinsp
end
twiny(plt::Plots.Plot=current()) = twiny(plt[1])

pct_formatter(yi) = "$(Int.(round.(yi * 100)))%"
blank_formatter(yi) = ""

function get_colormap()
    giallo = colorant"#F0BC42"
    rosso = colorant"#8E1F2F"
    grigio = colorant"#CACACC"
    nero = colorant"#000000"
    thirds = colorant"#31467A"
    orange = colorant"#F18101"
    return [rosso, giallo, thirds, grigio, orange, nero]
end