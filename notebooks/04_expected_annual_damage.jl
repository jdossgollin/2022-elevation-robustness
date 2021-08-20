### A Pluto.jl notebook ###
# v0.14.9

using Markdown
using InteractiveUtils

# ╔═╡ 82c386f3-c1a1-480e-b034-14034dc85ffe
begin
	using DrWatson
	using Plots
	using Unitful

	using DecisionProblem
end

# ╔═╡ 58422108-e041-11eb-10ce-13568cddf9bf
md"""
# Emulate Damages

This notebook gives a brief explanation of the model used to calculate expected annual damage as a function of ((house elevation) - (mean sea level)).
This model is already implemented in the `DecisionProblem` local package.
"""

# ╔═╡ f176c2b4-e441-4878-9753-b42ea6e4084b
md"""
The expected damages $d$ as a function of $h$, the difference between the elevation of the house and mean sea level, is the convolution of the probability density of storm surge $p(y')$ and the depth-damage function $d(y' - h)$:

$$\mathbb{E}(d | h) = \int p(y') \, d(y') \, dy'$$

This can be approximated as

$$\frac{1}{N} \sum_{i=1}^N d(y^*_i)$$

where $y^*_1, \ldots, y^*_N$ are draws from $p(y')$.
Approximating this integral many times over (eg, for each trajectory of sea level rise and each year) may be computationally costly.
Instead, we develop an emulator for this by calculating $\mathbb{E}(d | h)$ for many values of $h$ and then fitting a curve to this data.
"""

# ╔═╡ edf06404-68e7-476e-af48-9550f9cfd696
surge_mc = DecisionProblem.draw_surges(1_000_000);

# ╔═╡ caa24b9c-d816-4aaa-af3a-fb93d03a5897
depths = collect(-25:0.5:20)u"ft";

# ╔═╡ 411af8d9-362c-46af-8a01-42c89da9cc1d
xp_dmg = [
	# this function calculates the Monte Carlo integral
	DecisionProblem.expected_damage_frac(depth, surge_mc; key = :hazus)
	for depth in depths
];

# ╔═╡ 5f96d337-19df-4763-911b-577559fbf969
md"and calculate the expected (per year) structural damage for each"

# ╔═╡ b4fe8c62-c074-4369-aa86-36a7bf9b0564
begin
	depth_ft = ustrip.(u"ft", depths)
	hazus_emulator = get_expected_damage_emulator(:hazus)
	p = plot(
		depth_ft,
		100 * hazus_emulator.(depths),
		label="HAZUS",
		xlabel="House Elevation Relative to Mean Sea Level (ft)",
		ylabel="Expected Annual Damages, % House Value",
		marker=".",
		markersize=2,
		linewidth=2,
	)
	europa_emulator = get_expected_damage_emulator(:europa)
	plot!(
		p,
		depth_ft,
		100 * europa_emulator.(depths),
		label="Europa",
		marker=".",
		markersize=2,
		linewidth=2,
	)
	savefig(plotsdir("expected_annual_damage.pdf"))
	p
end

# ╔═╡ 21a22b3c-6cfe-450a-b99b-943f7003a909
md"""
This gives us a very clean relationship: given a house's elevation relative to the gauge and the mean sea level, we can now very quickly approximate expected (per year) damages.
This will make optimization much faster!
Since we've used a million Monte Carlo simulations of storm surge for each height difference, our results will be numerically stable.

This also works for the Europa depth-damage model.
"""

# ╔═╡ Cell order:
# ╟─58422108-e041-11eb-10ce-13568cddf9bf
# ╠═82c386f3-c1a1-480e-b034-14034dc85ffe
# ╟─f176c2b4-e441-4878-9753-b42ea6e4084b
# ╠═edf06404-68e7-476e-af48-9550f9cfd696
# ╠═caa24b9c-d816-4aaa-af3a-fb93d03a5897
# ╠═411af8d9-362c-46af-8a01-42c89da9cc1d
# ╟─5f96d337-19df-4763-911b-577559fbf969
# ╠═b4fe8c62-c074-4369-aa86-36a7bf9b0564
# ╟─21a22b3c-6cfe-450a-b99b-943f7003a909
