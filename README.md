# ReEDS2PRAS
ReEDS2PRAS

The purpose of ReEDS2PRAS is to translate completed ReEDS runs into PRAS systems ready for probabilistic resource adequacy analysis.

## Introduction
ReEDS2PRAS can be added in the Julia REPL from its git URL

`] add https://github.nrel.gov/PRAS/ReEDS2PRAS.git`

Assuming  the user can access the julia REPL and clone ReEDS2PRAS, the other necessary component for running standalone ReEDS2PRAS is a completed ReEDS run. More information about running ReEDS is available at `https://github.nrel.gov/ReEDS/ReEDS-2.0`.

## Basic Usage
If you have a completed ReEDS run and a REPL with ReEDS2PRAS (`using ReEDS2PRAS`), an example of running ReEDS2PRAS is provided below

```
using ReEDS2PRAS

# directory where user has ReEDS. This is used to find the input EIA database for disaggregating generators.
nems_path = joinpath("/projects/ntps/llavin/ReEDS-2.0")
# path to completed ReEDS runs
reeds_filepath = joinpath("/projects/ntps/llavin/ReEDS-2.0/runs/ntpsrerun_Xlim_DemHi_90by2035EarlyPhaseout__core")
solve_year = 2035
weather_year = 2012 # must be 2007-2013
timesteps = 8760

# returns a parameterized PRAS system
pras_system = ReEDS2PRAS.reeds_to_pras(reeds_filepath, solve_year, nems_path,
                                       timesteps, weather_year)
```

This will save out a pras system to the variable `pras_system` from the ReEDS2PRAS run. The user can also save a PRAS system to a specific location using `PRAS.savemodel(pras_system, joinpath("MYPATH"*".pras")`. The saved PRAS system may then be read in by other tools like PRAS Analytics (`https://github.nrel.gov/PRAS/PRAS-Analytics`) for further analysis, post-processing, and plotting.

### This can also be run through the command line:

```
julia --project=. src/ReEDS2PRAS.jl <reeds_filepath> <solve_year> <nems_path> <timesteps> <weather_year> <output_filepath>
```

The `output_filepath` specifies the location to save the pras model.


## Acknowledgements
The developers are Luke Lavin, Surya Dhulipala, and Brandon Benton. They acknowlege Gord Stephen for many helpful comments for improving this package.