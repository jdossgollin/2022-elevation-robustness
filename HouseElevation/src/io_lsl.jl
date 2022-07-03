using CSV
using DataFrames
using Interpolations
using NetCDF
using Unitful

struct LSLSim{I<:Integer,T<:AbstractFloat}
    years::UnitRange{I}
    lsl::Vector{<:Unitful.Length}
    rcp::T
    model::AbstractString
end

"""A convenient function to convert decadal to annual data"""
function decadal_to_annual(year_in, val_in, year_out)
    extrap = LinearInterpolation(year_in, val_in; extrapolation_bc=Line())
    return extrap.(year_out)
end

"""
Read in the subsidence data from Kopp et al (2014)

Kopp, R. E., Horton, R. M., Little, C. M., Mitrovica, J. X., Oppenheimer, M., Rasmussen, D. J., et al. (2014). Probabilistic 21st and 22nd century sea-level projections at a global network of tide-gauge sites. Earth’s Future, 2(8), 383-406. https://doi.org/10.1002/2014ef000239
"""
function get_K14_subsidence(;
    syear::Union{Int,Missing}=missing, eyear::Union{Int,Missing}=missing
)
    fname = data_dir("raw", "msl", "LSLProj_bkgd_299_rcp26.csv")
    subsidence = Array(DataFrames.DataFrame(CSV.File(fname)))
    years = collect(2010:10:2200)
    if ismissing(syear)
        syear = 2010
    end
    if ismissing(eyear)
        eyear = 2200
    end
    subsidence_annual =
        hcat(
            [
                decadal_to_annual(years, row, collect(syear:eyear)) for
                row in eachrow(subsidence)
            ]...,
        ) .* 1u"cm"
    return subsidence_annual
end

"""
These simulations use the BRICK model:

Wong, T. E., Bakker, A. M. R., Ruckert, K., Applegate, P., Slangen, A. B. A., & Keller, K. (2017). BRICK v0.2, a simple, accessible, and transparent model framework for climate and regional sea-level projections. Geoscientific Model Development, 10(7), 2741-2760. https://doi.org/10.5194/gmd-10-2741-2017

As noted in Ruckert et al

Ruckert, K. L., Srikrishnan, V., & Keller, K. (2019). Characterizing the deep uncertainties surrounding coastal flood hazard projections: A case study for Norfolk, VA. Scientific Reports, 9(1), 1-12. https://doi.org/10.1038/s41598-019-47587-6

these data do not have subsidence, so we add in the subsidence from Kopp et al

Kopp, R. E., Horton, R. M., Little, C. M., Mitrovica, J. X., Oppenheimer, M., Rasmussen, D. J., et al. (2014). Probabilistic 21st and 22nd century sea-level projections at a global network of tide-gauge sites. Earth’s Future, 2(8), 383-406. https://doi.org/10.1002/2014ef000239
"""
function get_lsl_brick(; syear::Int=2022, eyear::Int=2200)
    fnames = [
        "BRICK_NOfastDynamics_SP_20Nov2018.nc",
        "BRICK_SewellsPoint_FastDynamics_20Nov2018.nc",
    ]
    models = ["BRICK Slow", "BRICK Fast"]

    years = Int.(NetCDF.ncread(data_dir("raw", "msl", first(fnames)), "time_proj"))
    idx_syear = findfirst(years .== syear)

    if ismissing(eyear)
        eyear = maximum(years)
    end
    idx_eyear = findfirst(years .== eyear)

    # read in subsidence data
    subsidence = get_K14_subsidence(; syear=syear, eyear=eyear)
    n_subsidence = size(subsidence, 2)

    s = []
    for (fname_short, model) in zip(fnames, models)
        fname = data_dir("raw", "msl", fname_short)

        # read in the raw data and convert to meters
        lsl_26 = NetCDF.ncread(fname, "LocalSeaLevel_RCP26") .* 1u"m"
        lsl_45 = NetCDF.ncread(fname, "LocalSeaLevel_RCP45") .* 1u"m"
        lsl_60 = NetCDF.ncread(fname, "LocalSeaLevel_RCP60") .* 1u"m"
        lsl_85 = NetCDF.ncread(fname, "LocalSeaLevel_RCP85") .* 1u"m"

        # convert to LSLSim, noting that the underlying data is in meters
        scenarios = [
            (rcp=2.6, model=model, data=lsl_26),
            (rcp=4.5, model=model, data=lsl_45),
            (rcp=6.0, model=model, data=lsl_60),
            (rcp=8.5, model=model, data=lsl_85),
        ]
        si = vcat(
            [
                [
                    LSLSim(
                        syear:eyear,
                        let
                            brick_lsl = scenario.data[idx_syear:idx_eyear, i]
                            subsid = subsidence[:, rand(1:n_subsidence)]
                            net = brick_lsl .- subsid
                            net .- first(net)
                        end,
                        scenario.rcp,
                        scenario.model,
                    ) for i in 1:size(scenario.data, 2)
                ] for scenario in scenarios
            ]...,
        )
        # normalize the data so that it's relative to syear
        push!(s, si)
    end
    s = vcat(s...)
    return s
end

"""
Get probabilistic sea level rise data from Kopp et al (2017)

Kopp, R. E., DeConto, R. M., Bader, D. A., Hay, C. C., Horton, R. M., Kulp, S., et al. (2017). Evolving understanding of Antarctic ice-sheet physics and ambiguity in probabilistic sea-level projections. Earth’s Future, 5(12), 1217-1233. https://doi.org/10.1002/2017ef000663

This paper provides results using two ice sheet models: DP16 and K14.

DeConto, R. M., & Pollard, D. (2016). Contribution of Antarctica to past and future sea-level rise. Nature, 531(7596), 591-597. https://doi.org/10.1038/nature17145

Kopp, R. E., Horton, R. M., Little, C. M., Mitrovica, J. X., Oppenheimer, M., Rasmussen, D. J., et al. (2014). Probabilistic 21st and 22nd century sea-level projections at a global network of tide-gauge sites. Earth’s Future, 2(8), 383-406. https://doi.org/10.1002/2014ef000239
"""
function get_lsl_k17(; syear::Int=2022, eyear::Int=2200)
    scenarios = [
        (rcp=2.6, model="K14", fname="LSLproj_MC_299_rcp26.csv"),
        (rcp=4.5, model="K14", fname="LSLproj_MC_299_rcp45.csv"),
        (rcp=6.0, model="K14", fname="LSLproj_MC_299_rcp60.csv"),
        (rcp=8.5, model="K14", fname="LSLproj_MC_299_rcp85.csv"),
        (rcp=2.6, model="DP16", fname="LSLproj_MC_DP16_SEW_299_rcp26.csv"),
        (rcp=4.5, model="DP16", fname="LSLproj_MC_DP16_SEW_299_rcp45.csv"),
        (rcp=6.0, model="DP16", fname="LSLproj_MC_DP16_SEW_299_rcp60.csv"),
        (rcp=8.5, model="DP16", fname="LSLproj_MC_DP16_SEW_299_rcp85.csv"),
    ]
    years = collect(2010:10:2200)
    if ismissing(syear)
        syear = minimum(years)
    end
    if ismissing(eyear)
        eyear = maximum(years)
    end
    years_annual = syear:eyear

    s = []
    for scenario in scenarios
        msl = Array{Float64}(
            select(
                DataFrame(CSV.File(data_dir("raw", "msl", scenario.fname))), string.(years)
            ),
        )
        msl_annual =
            hcat(
                [
                    decadal_to_annual(years, row, collect(syear:eyear)) for
                    row in eachrow(msl)
                ]...,
            ) .* 1u"cm"
        si = [
            LSLSim(
                years_annual,
                let
                    msl_i = msl_annual[:, i]
                    msl_i .- first(msl_i)
                end,
                scenario.rcp,
                scenario.model,
            ) for i in 1:size(msl_annual, 2)
        ]
        push!(s, si)
    end
    return vcat(s...)
end

"""
Get all the LSL simulations
"""
function get_lsl(; syear::Int=2022, eyear::Int=2200)
    return vcat(
        get_lsl_brick(; syear=syear, eyear=eyear), get_lsl_k17(; syear=syear, eyear=eyear)
    )
end

function get_year_data(t::LSLSim, y::Int)
    idx = findfirst(t.years .== y)
    return t.lsl[idx]
end

function get_year_data(ts::Vector{<:LSLSim}, y::Int)
    idx = findfirst(first(ts).years .== y)
    return [t.lsl[idx] for t in ts]
end
