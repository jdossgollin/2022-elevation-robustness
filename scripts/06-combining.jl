using LaTeXStrings

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
        xlabel="Data (x)",
        legend=:right,
        bottom_margin=5Plots.mm,
        left_margin=5Plots.mm,
        size=(800, 400),
    )
    scatter!(p, ψ, zeros(size(ψ)); label="Samples", color=colors[2], markersize=5)

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
            color=colors[3],
            linestyle=:dash,
        )
        plot!(p, [xmax, xmax], [0, ymax]; label=false, color=colors[3], linestyle=:dash)
        if i < length(ψ)
            plot!(
                p2, [xmax, xmax], [0, ymax]; color=colors[2], label=false, linestyle=:dash
            )
        end
        plot!(
            p,
            [x_arrow[i], x_arrow[i]],
            [ymax_prev, ymax];
            color=colors[4],
            linewidth=3,
            label=false,
        )
        annotate!(
            p,
            [x_arrow[i] + 0.2],
            [(ymax_prev + ymax) / 2],
            text(L"$w_%$i$", :center, 10; color=colors[4]),
        )
        ymax_prev = ymax
    end
    savefig(p, plots_dir("grid-sketch.pdf"))
    return p
end
