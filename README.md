# 2021-elevation-robustness

This code is developed by James Doss-Gollin.
If you use this code, please cite our paper:

```bibtex
TBD
```

## Setup

### Install Julia

All code has been developed using Julia 1.6, specifically the long term support (LTS) version 1.6.5.
You can download this version of Julia at [https://julialang.org/downloads/](https://julialang.org/downloads/).

### Install Packages

This repository uses a Julia environment.
To download all required packages, run the following in a Julia session

```julia
using Pkg
Pkg.instantiate()
```

### Running

In a Julia session (assuming you have navigated to this project directory) run

```julia
include("scripts/main.jl")
```

to run everything.
You can also run the lines in `main.jl` interactively.
Alternatively, you can run all analysis as a script from the command line with

```bash
julia --project scripts/main.jl
```
