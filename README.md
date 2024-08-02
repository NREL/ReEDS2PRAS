# ReEDS2PRAS
## Introduction

The purpose of ReEDS2PRAS is to translate completed ReEDS runs into PRAS systems ready for probabilistic resource adequacy analysis.

### Installation

ReEDS2PRAS can be installed in the Julia REPL from the `NRELInternal` package
registry. The (public) `NREL` package registry is also required as that's
where PRAS is listed.

If one or both of those registries haven't already been added, enter the
package management prompt (hit `]` from the standard Julia prompt) and run:

```
registry update
registry add https://github.com/NREL/JuliaRegistry.git
registry add https://github.nrel.gov/JuliaLang/InternalRegistry.git
```

The first command just makes sure you also add the default `General` Julia
registry. The second and third add the public and private NREL registries,
respectively. Once these registries are set up on your machine you don't need
to repeat this process  - you can always check which registries you have set up
by running `registry status`.

To install the actual package, run `add ReEDS2PRAS` while still in the
package management prompt. Specific versions of the package can also be
installed, e.g.: `add ReEDS2PRAS@v0.2.5` (see the compatibility table below).

Once this is finished, hit backspace to return to the standard Julia prompt.

Assuming  the user can access the julia REPL and install ReEDS2PRAS, the other
necessary component for running standalone ReEDS2PRAS is a completed ReEDS run.
More information about running ReEDS is available at
`https://github.nrel.gov/ReEDS/ReEDS-2.0`.

### Compatibility Table

This won't cover every possible version combination, but should at least
provide a starting point for determining which version of the package to use:

 For this version of ReEDS... | ...use this version of ReEDS2PRAS
------------------------------|----------------------------------
2023.0.0 | 0.2.5
2023.40.3 | 0.3.0

Unless indicated otherwise, ReEDS2PRAS should be expected to produce
a PRAS v0.6 system.

## Basic Usage
If you have a completed ReEDS run and a REPL with ReEDS2PRAS (`using ReEDS2PRAS`), an example of running ReEDS2PRAS is provided below

```
using ReEDS2PRAS

reedscase = "/projects/ntps/llavin/ReEDS-2.0/runs/ntpsrerun_Xlim_DemHi_90by2035EarlyPhaseout__core" # path to completed ReEDS run
solve_year = 2035 #need ReEDS Augur data for the input solve year
weather_year = 2012 # must be 2007-2013
timesteps = 8760
user_descriptors = "your_user_descriptors_json_location_here" # Optional - if not passed uses default values

# returns a parameterized PRAS system
pras_system = ReEDS2PRAS.reeds_to_pras(reedscase, solve_year, timesteps, weather_year, user_descriptors = user_descriptors)
```

This will save out a pras system to the variable `pras_system` from the ReEDS2PRAS run. The user can also save a PRAS system to a specific location using `PRAS.savemodel(pras_system, joinpath("MYPATH"*".pras")`. The saved PRAS system may then be read in by other tools like PRAS Analytics (`https://github.nrel.gov/PRAS/PRAS-Analytics`) for further analysis, post-processing, and plotting.

## Multi-year usage

ReEDS2PRAS can be run for multiple weather years of a completed ReEDS run by passing more than 8760 hourly timestamps. For example running all 7 weather years can be accomplished as in the below example

```
using ReEDS2PRAS

# path to completed ReEDS run
reedscase = "/projects/ntps/llavin/ReEDS-2.0/runs/ntpsrerun_Xlim_DemHi_90by2035EarlyPhaseout__core"
solve_year = 2035 #need ReEDS Augur data for the input solve year
weather_year = 2007 # must be 2007-2013
timesteps = 61320
user_descriptors = "your_user_descriptors_json_location_here" # Optional - if not passed uses default values

# returns a parameterized PRAS system
pras_system = ReEDS2PRAS.reeds_to_pras(reedscase, solve_year, timesteps, weather_year, user_descriptors = user_descriptors)
```

Importantly, the timesteps count from the first hour of the first `weather_year`, so the user must input `2007` as the `weather_year` to run all 61320 hourly timesteps.

### This can also be run through the command line:

```
julia --project=. bin/run.jl <reedscase> <solve_year> <timesteps> <weather_year> <output_filepath>
```

The `output_filepath` specifies the location to save the pras model.

## Contributing Guidelines

It is always a good practice to follow the Julia Style Guide (`https://docs.julialang.org/en/v1/manual/style-guide/`). 
Please make sure you format your code to follow our guidelines using the snippet below before you open a PR:
```
julia  -e 'using Pkg; Pkg.add("JuliaFormatter"); using JuliaFormatter; include(".github/workflows/formatter-code.jl")'
```
**NOTE: You have to run the snippet above at the repo folder level.

## Acknowledgements
The developers are Luke Lavin, Surya Dhulipala, Brandon Benton, Patrick Brown, and Hari Sundar. They acknowledge Gord Stephen for many helpful comments for improving this package.
