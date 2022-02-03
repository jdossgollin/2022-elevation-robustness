using Distributions

function add_panel_letters!(panels; fontsize=9)
    letters = 'A':'Z'
    for (p, letter) in zip(panels, letters)
        annotate!(p, (0.05, 0.95), text("($letter)", :left, fontsize))
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
