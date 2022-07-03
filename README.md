# Repository for "A subjective Bayesian framework for synthesizing deep uncertainties in climate risk management"

This is the GitHub repository for the paper "A subjective Bayesian framework for synthesizing deep uncertainties in climate risk management" by James Doss-Gollin and Klaus Keller.
This code is developed by James Doss-Gollin.
If you use this code, please cite our paper:

```bibtex
TBD
```

## Reproducibility

This section provides guidance on reproducing our results.

We do not purport to achieve bitwise reproducibility, so your results may vary slightly, but you should be able to follow the above steps to get highly similar results.
If you are unable to reproduce our results, please open an Issue using the issues tab above.
The more detail you are able to provide, the more readily we can track down any potential problems.

### Install Julia

All code has been developed using Julia 1.6, specifically the long term support (LTS) version 1.6.5.
You can download this version of Julia at [https://julialang.org/downloads/](https://julialang.org/downloads/).

### Install Packages

This repository uses a Julia environment.
To download all required packages, run the following in a Julia session

```julia-repl
using Pkg
Pkg.instantiate()
```

*or* open `Pkg` mode in the REPL (by typing `]`) and run `instantiate`.

### Running

In a Julia REPL (assuming you have navigated to this project directory) run

```julia-repl
include("scripts/main.jl")
```

to run everything.
You can also run the lines in `main.jl` interactively.
Alternatively, you can run all analysis as a script from the command line with

```bash
julia --project scripts/main.jl
```

## Code organization

This section provides some additional details on how code is organized in this repository.

All code is written in Julia and can be run through one file.
The `HouseElevation` folder provides a local module.
Essentially, this module contains abstract code for a generic house elevation problem.

The results of this paper are produced by running the files in `scripts/` in numerical order.
The `main.jl` script accomplishes this, along with all necessary imports.
As suggested above, the best way to run our codes is to run this `main.jl` file.
The `scripts/` directory also includes `plotutils.jl`, which provides some functions used in the other scripts.
We did not feel that creating a package for these functions added value.

The folder `data/raw` contains required input data that we provide.
As you run scripts, intermediate results will be cached in `data/external/` or `data/processed/`.
If you delete these intermediate files, they will be re-created as needed.

The folder `papers/` contains LaTeX code used to create the manuscript of this file, as well as a poster presented at AGU 2021.

The folder `plots/` contains the figures and tables used in our paper.
If you re-run the codes, these may change slightly.
You are welcome to reuse these figures, but please cite our work and note that copyright of these figures technically belongs to the journal.

The folder `tikz/` contains some code to produce figures directly in LaTeX.

## DOI Releases

This repository contains a live repository of our code, incorporating updates and suggestions made over time.
We have used Zenodo to archive the precise versions of our code used to generate journal submissions.
Please see:

1. (Empty for now)
