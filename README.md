# A subjective Bayesian framework for synthesizing deep uncertainties in climate risk management

[![DOI](https://zenodo.org/badge/357754608.svg)](https://zenodo.org/badge/latestdoi/357754608)

This is the GitHub repository for the paper _**A subjective Bayesian framework for synthesizing deep uncertainties in climate risk management**_ by James Doss-Gollin and Klaus Keller.
This code is developed by James Doss-Gollin.
If you use our code, please cite our paper:

```bibtex
TBD
```

## Reproducibility

This section provides guidance on reproducing our results.
We do not purport to achieve bitwise reproducibility, so your results may vary slightly, but you should be able to follow the above steps to get highly similar results.
If you are unable to reproduce our results, please [open an Issue](https://github.com/jdossgollin/2022-elevation-robustness/issues).
You may also send a Tweet to [@jdossgollin](https://twitter.com/jdossgollin).
The more detail you are able to provide, the more readily we can track down any potential problems.

### Install Julia

All code has been developed using Julia 1.6, specifically the long term support (LTS) version 1.6.5.
You can download this version of Julia at [https://julialang.org/downloads/](https://julialang.org/downloads/).

### Install Packages

This repository uses a Julia environment.
To download all required packages, run the following in a Julia session (i.e, the REPL).

1. **Open a Julia REPL**. If you've never used Julia before, the [course webpage](https://jdossgollin.github.io/environmental-data-science/) has some helpful resources.
1. **Open `Pkg` mode** in the REPL (by typing `]`) and make sure you see something that looks like

    ```julia-repl
    (2022-elevation-robustness) pkg>
    ```

1. **Activate** the environment. Now that you are in `Pkg` mode, run `activate .` (and hit Enter)
1. **Instantiate** to set up the packages. To do this, just run `instantiate` (and hit Enter)
1. **Leave `Pkg` mode** by hitting Backspace

If you don't want to use `Pkg` mode, you can run the following in the Julia REPL

```julia-repl
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

### Running

In a Julia REPL, make sure that you have activated the project environment (see above) and run

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

1. [10.5281/zenodo.6799457](doi.org/10.5281/zenodo.6799457): code for generating preprint
