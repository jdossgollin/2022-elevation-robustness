# 2021-elevation-robustness

This code is developed by James Doss-Gollin.
If you use this code, please cite our paper:

```bibtex
TBD
```

## Setup

### Install Julia

YOU MUST USE JULIA 1.6.
We recommend the long term support (LTS) version 1.6.5.
If you try to install Julia 1.7 or above, you will run into errors!

### Install Packages

This repository uses a Julia environment.
To download all required packages, run the following in a Julia session

```julia
using Pkg
Pkg.instantiate()
```

### Running

In a Julia session (assuming you are in this project directory), you can run

```julia
include("scripts/main.jl")
```

to run everything.
You can also run the lines in `main.jl` interactively.
Alternatively, you can run all analysis as a script from the command line

```bash
julia --project scripts/main.jl
```

## Abstract

Mainstreaming climate adaptation requires subjective syntheses of deep uncertainties

Keller & Doss-Gollin

The policy and engineering tools used to inform climate adaptation etc rely on standards and cost-benefit analyses that depend on risk probabilities.
One example is house elevation: the recommended elevation is a standard (BFE plus a foot) and funding decisions rely on CBA.
Under stationarity, the probabilities for this can almost be agreed upon.
Under nonstationarity, there are deep uncertainties.
We show that the way these deep uncertainties are synthesized into CBA matters a lot.
Thus, static decisions are subject to deep uncertainties
