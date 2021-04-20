### A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

# ╔═╡ cfaa1f89-3fa7-4ca3-b32b-88c68d2e2081
using DrWatson

# ╔═╡ 04263916-8f87-4614-bb5d-a0ef1a1d7fe9
using PlutoUI; PlutoUI.TableOfContents(title="Outline")

# ╔═╡ b1709ffe-44c7-4f5c-a254-4a921717163e
using Plots, Unitful, ColorSchemes, UnitfulRecipes

# ╔═╡ d255b46c-1c29-42ff-a94e-db570fb029e8
using NorfolkFloods, HouseElevation

# ╔═╡ fcfb26dd-1767-49f7-ae97-259e09790be5
md"# Raw Data Plots

In this notebook we will plot some constants and fixed quantities from our models

## Front Matter

Start with package imports and definitions of constants
"

# ╔═╡ 386000ba-3b50-4f55-b674-065149a2eddf
@quickactivate "2021-elevation-robustness"

# ╔═╡ 68e96bb1-f5ef-47e3-b293-d3c915a2e558
colors = ColorSchemes.tab10;

# ╔═╡ d03aeb0d-ab72-49f3-8352-c5ddc4d475fb
md"## Historic Storms

Let's start by plotting historic annual  maximum storm surges.
We will annotate these with some information (from Wikipedia, NWS, etc) on what caused these surges.

We plot surge plus sea level.
"

# ╔═╡ 7acfee42-a1ff-11eb-2ddb-93669963b4a6
function plot_historic_annmax_flood()

    annual = get_norfolk_annual()
    surge_ft = (annual.msl + annual.surge) .|> x -> ustrip(u"ft", x)

    p = plot(
        annual.year,
        surge_ft,
        label = "",
        ylabel = "Annual Maximum Storm Surge (ft)",
        marker = :circle,
        size = (500, 300),
        title = "NOAA $(annual.gage_id): $(annual.gage_name)",
        titlefontsize = 10,
        labelfontsize = 9,
    )

    for (year, name, Δx, Δy, is_tc) in zip(
        [1933, 2003, 2015, 2009, 1962],
        ["Chesapeake–Potomac", "Isabel", "Irene", "Nor'Ida", "Ash Wednesday"],
        [12.5, -7.5, 2.5, 0, 10],
        [0.75, 0.5, 0.5, 1, 0.75],
        [true, true, true, false, false],
    )
        yobs = ustrip(surge_ft[findfirst(annual.year .== year)])
        x0 = year + Δx
        y0 = yobs + Δy
        color = is_tc ? colors[2] : colors[1]
        plot!(
            [x0, year],
            [y0, yobs],
            arrow = :closed,
            color = color,
            linewidth = 1,
            label = "",
        )
        scatter!(
            p,
            [x0],
            [y0],
            markercolor = :white,
            label = false,
            markerstrokecolor = :white,
            markersize = 10,
        )
        annotate!(p, x0, y0, text(name, :center, 7, color = color))
    end
    return p
end

# ╔═╡ 972e6a82-819d-4d1e-b240-5a8e88fdcb61
annmax_flood = plot_historic_annmax_flood()

# ╔═╡ 6d1ef7b4-653e-4046-a49c-884c152c5489
savefig(annmax_flood, plotsdir("ts_historic_anmax_flood.pdf"))

# ╔═╡ 3e0aafd4-5777-4fbe-ae87-17d28f1bd91b
obs = get_norfolk_annual();

# ╔═╡ 6cd4d7b1-eae0-4a17-92c0-638135c994c9
md"""
The largest flood occurred in 1933 as part of the Chesapeake–Potomac hurricane.
Since then there has been significant sea level rise; if that event occurred with current sea level it might be much higher
Specifically, if it occurred with the maximum observed mean sea level then it would be

$(uconvert(u"ft", maximum(obs.surge) + maximum(obs.msl)))
"""

# ╔═╡ 785eb441-acb7-46dd-9ce3-61847a789bcf
md"""
## Elevation Cost

Now let's visualize the cost of elevating a house -- this is another fixed quantity that we can examine
"""

# ╔═╡ a0581141-df6b-422f-9917-14fdb4fa05ca
md"""
The most striking thing about the elevation cost is that it is very gentle for small increases.
Thus, if it makes sense to elevate a little bit, it probably makes sense to elevate a lot.

This curve is very important, but we are relying on rather shaky data (see Zarekarizi et al, 2020); using more accurate elevation cost data would *definitely change results*

- Zarekarizi, M., Srikrishnan, V., & Keller, K. (2020). Neglecting uncertainties biases house-elevation decisions to manage riverine flood risks. Nature Communications, 11(1), 5361. https://doi.org/10.1038/s41467-020-19188-9
"""

# ╔═╡ 20cb9bcb-e57c-46fa-9c43-75685b85eef7
savefig(elevation_cost, plotsdir("elevation_cost.pdf"))

# ╔═╡ 436745c1-6865-4d67-9604-9452dd05870a
md"""
## Depth Damage

Now we plot the depth-damage functions that we are using
"""

# ╔═╡ 6929dab9-204a-41cc-91b4-391ab33c0fea
function plot_depth_damage()

    depths = (-2:0.1:30)u"ft"
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

- Rözer, V., Kreibich, H., Schröter, K., Müller, M., Sairam, N., Doss-Gollin, J., et al. (2019). Probabilistic models significantly reduce uncertainty in Hurricane Harvey pluvial flood loss estimates. Earth’s Future, 7(4). https://doi.org/10.1029/2018ef001074
- Rözer, V., Müller, M., Bubeck, P., Kienzler, S., Thieken, A., Pech, I., et al. (2016). Coping with pluvial floods by private households. Water, 8(7), 304. https://doi.org/10.3390/w8070304
- Johnson, J. F. (2000). Generic depth-damage relationships (Economic Guidance Memorandum No. EGM 01-03). United States Army Corps of Engineers Planning and Policy Division. Retrieved from https://planning.erdc.dren.mil/toolbox/library/EGMs/egm01-03.pdf
- de Moel, H., van Vliet, M., & Aerts, J. C. J. H. (2014). Evaluating the effect of flood damage-reducing measures: a case study of the unembanked area of Rotterdam, the Netherlands. Regional Environmental Change, 14(3), 895–908. https://doi.org/10.1007/s10113-013-0420-z
"

# ╔═╡ f3bfa252-9fea-4371-be86-94cafd166956
savefig(depth_damage, plotsdir("depth_damage.pdf"))

# ╔═╡ 23d3bad3-7b4d-46f4-91e7-1978411856ea
elevation_cost = plot_elevation_cost()

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

# ╔═╡ Cell order:
# ╠═fcfb26dd-1767-49f7-ae97-259e09790be5
# ╠═cfaa1f89-3fa7-4ca3-b32b-88c68d2e2081
# ╠═386000ba-3b50-4f55-b674-065149a2eddf
# ╠═04263916-8f87-4614-bb5d-a0ef1a1d7fe9
# ╠═b1709ffe-44c7-4f5c-a254-4a921717163e
# ╠═d255b46c-1c29-42ff-a94e-db570fb029e8
# ╠═68e96bb1-f5ef-47e3-b293-d3c915a2e558
# ╟─d03aeb0d-ab72-49f3-8352-c5ddc4d475fb
# ╠═7acfee42-a1ff-11eb-2ddb-93669963b4a6
# ╠═972e6a82-819d-4d1e-b240-5a8e88fdcb61
# ╠═6d1ef7b4-653e-4046-a49c-884c152c5489
# ╠═3e0aafd4-5777-4fbe-ae87-17d28f1bd91b
# ╟─6cd4d7b1-eae0-4a17-92c0-638135c994c9
# ╟─785eb441-acb7-46dd-9ce3-61847a789bcf
# ╠═b02ce2d7-2128-4d1d-93ee-0781ea61fbb9
# ╠═23d3bad3-7b4d-46f4-91e7-1978411856ea
# ╟─a0581141-df6b-422f-9917-14fdb4fa05ca
# ╠═20cb9bcb-e57c-46fa-9c43-75685b85eef7
# ╠═436745c1-6865-4d67-9604-9452dd05870a
# ╠═6929dab9-204a-41cc-91b4-391ab33c0fea
# ╠═dd50959a-cbca-4241-ac7d-28b797521df8
# ╟─9da6e034-f071-4445-82be-5fd6835aee43
# ╠═f3bfa252-9fea-4371-be86-94cafd166956
