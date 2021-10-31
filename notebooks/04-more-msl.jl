### A Pluto.jl notebook ###
# v0.16.4

using Markdown
using InteractiveUtils

# ╔═╡ 62956480-37bc-11ec-1901-13f49754cd1d
begin
	using Pkg
	Pkg.activate("..")
	using CSV
	using DataFrames
	using Downloads
	using DrWatson
end

# ╔═╡ 1db7b953-25b2-4969-8b1f-d4af37b9223c
md"Some more scenarios"

# ╔═╡ 5e18bac7-5217-4e42-b1fa-0630449173d9
function get_file(fname::String, url::String)
    if !isfile(fname)
		Downloads.download(url, fname)
	end
	return 1
end;

# ╔═╡ efa86325-f600-484c-abac-53ecc4bded59
sources = [
	(
		url = "https://github.com/scrim-network/local-coastal-flood-risk/raw/master/Data/LSLproj_MC_299_rcp26.csv",
		fname=datadir("processed", "MC_200_rcp26.csv")
	)
]

# ╔═╡ d54110c9-975e-4192-bab0-e44dc7baf6be
[get_file(source.fname, source.url) for source in sources]

# ╔═╡ b45d4406-4b4b-49c4-a77f-8d2ab5c05d7f
df = CSV.File(first(sources).fname) |> DataFrame;

# ╔═╡ 4d18f67e-0376-4a89-9013-d9a035e8b882
df[!, :sample] = 1:nrow(df)

# ╔═╡ c34e10df-739b-49b9-9ca4-caaea5e632d6
df

# ╔═╡ Cell order:
# ╠═62956480-37bc-11ec-1901-13f49754cd1d
# ╠═1db7b953-25b2-4969-8b1f-d4af37b9223c
# ╠═5e18bac7-5217-4e42-b1fa-0630449173d9
# ╠═efa86325-f600-484c-abac-53ecc4bded59
# ╠═d54110c9-975e-4192-bab0-e44dc7baf6be
# ╠═b45d4406-4b4b-49c4-a77f-8d2ab5c05d7f
# ╠═4d18f67e-0376-4a89-9013-d9a035e8b882
# ╠═c34e10df-739b-49b9-9ca4-caaea5e632d6
