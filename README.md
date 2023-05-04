# ReEDS2PRAS
ReEDS2PRAS

The purpose of ReEDS2PRAS is to translate completed ReEDS runs into PRAS systems ready for probabilistic resource adequacy analysis.

## Introduction
ReEDS2PRAS can be added in the Julia REPL from its git URL

`] add https://github.nrel.gov/PRAS/ReEDS2PRAS.git`

Assuming  the user can access the julia REPL and clone ReEDS2PRAS, the other necessary component for running standalone ReEDS2PRAS is a completed ReEDS run. More information about running ReEDS is available at `https://github.nrel.gov/ReEDS/ReEDS-2.0`.

### Compatibility Table

This won't cover every possible version combination, but should at least
provide a starting point for determining which version of the package to use:

 For this version of ReEDS... | ...use this version of ReEDS2PRAS
------------------------------|----------------------------------
2023.0.0 | 0.2.5

Unless indicated otherwise, ReEDS2PRAS should be expected to produce
a PRAS v0.6 system.

## Basic Usage
If you have a completed ReEDS run and a REPL with ReEDS2PRAS (`using ReEDS2PRAS`), an example of running ReEDS2PRAS is provided below

```
using ReEDS2PRAS

# path to completed ReEDS run
reedscase = joinpath("/projects/ntps/llavin/ReEDS-2.0/runs/ntpsrerun_Xlim_DemHi_90by2035EarlyPhaseout__core")
solve_year = 2035 #need ReEDS Augur data for the input solve year
weather_year = 2012 # must be 2007-2013
timesteps = 8760

# returns a parameterized PRAS system
pras_system = ReEDS2PRAS.reeds_to_pras(reedscase, solve_year, timesteps, weather_year)
```

This will save out a pras system to the variable `pras_system` from the ReEDS2PRAS run. The user can also save a PRAS system to a specific location using `PRAS.savemodel(pras_system, joinpath("MYPATH"*".pras")`. The saved PRAS system may then be read in by other tools like PRAS Analytics (`https://github.nrel.gov/PRAS/PRAS-Analytics`) for further analysis, post-processing, and plotting.

## Multi-year usage

ReEDS2PRAS can be run for multiple weather years of ra completed ReEDS run by passing more than 8760 hourly timestamps. For example running all 7 weather years can be accomplished as in the below example

```
using ReEDS2PRAS

# path to completed ReEDS run
reedscase = joinpath("/projects/ntps/llavin/ReEDS-2.0/runs/ntpsrerun_Xlim_DemHi_90by2035EarlyPhaseout__core")
solve_year = 2035 #need ReEDS Augur data for the input solve year
weather_year = 2007 # must be 2007-2013
timesteps = 61320

# returns a parameterized PRAS system
pras_system = ReEDS2PRAS.reeds_to_pras(reedscase, solve_year, timesteps, weather_year)
```

Importantly, the timesteps count from the first hour of the first `weather_year`, so the user must input `2007` as the `weather_year` to run all 61320 hourly timesteps.

### This can also be run through the command line:

```
julia --project=. bin/run.jl <reedscase> <solve_year> <reedspath> <timesteps> <weather_year> <output_filepath>
```

The `output_filepath` specifies the location to save the pras model.


## Acknowledgements
The developers are Luke Lavin, Surya Dhulipala, and Brandon Benton. They acknowlege Gord Stephen for many helpful comments for improving this package.
