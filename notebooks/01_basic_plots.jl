### A Pluto.jl notebook ###
# v0.14.9

using Markdown
using InteractiveUtils

# ╔═╡ 0c531ffe-381a-448a-90d2-ca1bd42d28c6
begin
    using DrWatson
	using PlutoUI
    using Plots
    using Unitful
    using ColorSchemes
    using UnitfulRecipes

    using NorfolkFloods
    using HouseElevation

    PlutoUI.TableOfContents()
end


# ╔═╡ fcfb26dd-1767-49f7-ae97-259e09790be5
md"# Raw Data Plots

In this notebook we will plot some constants and fixed quantities from our models

## Front Matter

Start with package imports and definitions of constants
"

# ╔═╡ 0422680b-afb4-4e5e-a1fe-421562b1d66f
colors = ColorSchemes.okabe_ito;

# ╔═╡ 785eb441-acb7-46dd-9ce3-61847a789bcf
md"""
## Elevation Cost

Now let's visualize the cost of elevating a house -- this is another fixed quantity that we can examine.
This plot should reproduce supplemental figure S3 of Zarekarizi et al (2020).

> Zarekarizi, M., Srikrishnan, V., & Keller, K. (2020). Neglecting uncertainties biases house-elevation decisions to manage riverine flood risks. Nature Communications, 11(1), 5361. https://doi.org/10.1038/s41467-020-19188-9
"""

# ╔═╡ b02ce2d7-2128-4d1d-93ee-0781ea61fbb9
function plot_elevation_cost()

    areas = [3_000, 2_000, 1_000]u"ft^2"
    elevations = (0:0.1:14)u"ft"

    p = plot(
        legend = :bottomright,
        xlabel = "Height Increase",
        ylabel = "Elevation Cost (x1000 \$USD)",
        ymin = 0,
        palette = :tab10,
        xticks = 0:2:14,
        size = (600, 400),
    )
    for area in areas
        house = HouseStructure(area, 185_000, 5u"ft") # last 2 don't matter here
        plot!(
            p,
            elevations,
            [elevation_cost(house, Δh) for Δh in elevations] ./ 1_000,
            label = "$area house",
            linewidth = 3,
        )
    end
    return p
end;

# ╔═╡ cf822e50-9eb3-42a9-ba6b-c2c9d8748c2e
elevation_cost_fig = plot_elevation_cost()

# ╔═╡ a0581141-df6b-422f-9917-14fdb4fa05ca
md"""
The most striking thing about the elevation cost is that it is very gentle for small increases.
Thus, if it makes sense to elevate a little bit, it probably makes sense to elevate a lot.

This curve is very important, but we are relying on rather shaky data (see Zarekarizi et al, 2020); using more accurate elevation cost data would *very likely change results*.
"""

# ╔═╡ 20cb9bcb-e57c-46fa-9c43-75685b85eef7
savefig(elevation_cost_fig, plotsdir("elevation_cost.pdf"))

# ╔═╡ 436745c1-6865-4d67-9604-9452dd05870a
md"""
## Depth Damage

Now we plot the depth-damage functions that we are using
"""

# ╔═╡ 6929dab9-204a-41cc-91b4-391ab33c0fea
function plot_depth_damage()

    depths = (-3:0.1:30)u"ft"
    house = HouseStructure(1000u"ft^2", 185_000, 0u"ft")

    p = plot(
        legend = :bottomright,
        xlabel = "Flood Depth Relative to House",
        ylabel = "Damages (Proportion of House Value)",
        ymin = 0,
        palette = :tab10,
        size = (700, 500),
    )
    for model in [:hazus, :europa]
        plot!(
            p,
            depths,
            [depth_damage_frac(house, depth, model) for depth in depths],
            label = string(model),
            linewidth = 3,
        )
    end
    return p
end;

# ╔═╡ dd50959a-cbca-4241-ac7d-28b797521df8
depth_damage = plot_depth_damage()

# ╔═╡ 9da6e034-f071-4445-82be-5fd6835aee43
md"
There are large disparities between the HAZUS and EUROPA datasets, which is unsurprising (see below references for some discussion of data used.
It would certainly be advisable to use depth-damage functions better tailored to a particular house.
*Engineers should be aware that this design choice has a big effect and is very hard to measure*.

> Rözer, V., Kreibich, H., Schröter, K., Müller, M., Sairam, N., Doss-Gollin, J., et al. (2019). Probabilistic models significantly reduce uncertainty in Hurricane Harvey pluvial flood loss estimates. Earth’s Future, 7(4). https://doi.org/10.1029/2018ef001074

> Rözer, V., Müller, M., Bubeck, P., Kienzler, S., Thieken, A., Pech, I., et al. (2016). Coping with pluvial floods by private households. Water, 8(7), 304. https://doi.org/10.3390/w8070304

> Johnson, J. F. (2000). Generic depth-damage relationships (Economic Guidance Memorandum No. EGM 01-03). United States Army Corps of Engineers Planning and Policy Division. Retrieved from https://planning.erdc.dren.mil/toolbox/library/EGMs/egm01-03.pdf

> de Moel, H., van Vliet, M., & Aerts, J. C. J. H. (2014). Evaluating the effect of flood damage-reducing measures: a case study of the unembanked area of Rotterdam, the Netherlands. Regional Environmental Change, 14(3), 895–908. https://doi.org/10.1007/s10113-013-0420-z
"

# ╔═╡ f3bfa252-9fea-4371-be86-94cafd166956
savefig(depth_damage, plotsdir("depth_damage.pdf"))

# ╔═╡ Cell order:
# ╟─fcfb26dd-1767-49f7-ae97-259e09790be5
# ╠═0c531ffe-381a-448a-90d2-ca1bd42d28c6
# ╠═0422680b-afb4-4e5e-a1fe-421562b1d66f
# ╟─785eb441-acb7-46dd-9ce3-61847a789bcf
# ╠═b02ce2d7-2128-4d1d-93ee-0781ea61fbb9
# ╠═cf822e50-9eb3-42a9-ba6b-c2c9d8748c2e
# ╟─a0581141-df6b-422f-9917-14fdb4fa05ca
# ╠═20cb9bcb-e57c-46fa-9c43-75685b85eef7
# ╟─436745c1-6865-4d67-9604-9452dd05870a
# ╠═6929dab9-204a-41cc-91b4-391ab33c0fea
# ╠═dd50959a-cbca-4241-ac7d-28b797521df8
# ╟─9da6e034-f071-4445-82be-5fd6835aee43
# ╠═f3bfa252-9fea-4371-be86-94cafd166956
