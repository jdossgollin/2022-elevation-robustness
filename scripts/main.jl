using Revise
using DrWatson
@quickactivate "2021-elevation-robustness"

using PlutoSliderServer
using Pluto

function clean_outputs()
    for dir in [datadir("processed"), plotsdir(), projectdir("notebooks")]
        if isdir(dir)
            rm(dir; recursive = true)
        end
        mkdir(dir)
    end
end

function get_all_scripts()
    return scriptsdir.([
        "01_basic_plots.jl",
        "02_surge_model.jl",
        "03_plot_historical.jl",
        "04_expected_annual_damage.jl",
        "05_probabilistic_msl.jl",
    ])
end

function run_pluto()
    Pluto.run()
end

function convert_notebooks()
    notebooks = get_all_scripts()

    PlutoSliderServer.export_notebook.(
        notebooks;
        Export_output_dir = projectdir("notebooks"),
    )
end

function source_scripts()
    scripts = get_all_scripts()
    for script in scripts
        @info "Now in $script"
        include(script)
    end
end

# run one of the commented commands
# source_scripts()
# run_pluto()
# convert_notebooks()