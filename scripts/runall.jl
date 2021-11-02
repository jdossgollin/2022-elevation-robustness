using DrWatson
using PlutoSliderServer: export_notebook
using Pluto

NB_DIR = projectdir("notebooks")
OUTPUT_DIR = projectdir("nb_out")


"""Helper function to get the extension of a file"""
function file_extension(filename)
    return filename[findlast(isequal('.'), filename):end]
end

"""Clean up all processed data, plots, and outputs"""
function clean_outputs()
    for dir in [datadir("processed"), plotsdir(), OUTPUT_DIR]
        if isdir(dir)
            rm(dir; recursive = true)
        end
        mkdir(dir)
    end
end

"""Helper function to get all the notebooks"""
function get_all_notebooks()
    return [
        joinpath(NB_DIR, file) for file in readdir(NB_DIR) if file_extension(file) == ".jl"
    ]
end

"""run Pluto"""
function run_pluto(; port = 8888)
    Pluto.run(port = port)
end

"""run all notebooks and convert to HTML files"""
function convert_notebooks()
    if !isdir(OUTPUT_DIR)
        mkpath(OUTPUT_DIR)
    end
    for nb in get_all_notebooks()
        export_notebook(nb, Export_output_dir = OUTPUT_DIR)
    end
end

# run one of the following
# clean_outputs()
# run_pluto()
# convert_notebooks()
