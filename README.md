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

NP = joinpath("/projects/ntps/llavin/ReEDS-2.0") #directory where user has ReEDS. This is used to find the input EIA database for disaggregating generators.
runpath = joinpath("/projects/ntps/llavin/ReEDS-2.0/runs/ntpsrerun_Xlim_DemHi_90by2035EarlyPhaseout__core") #path to completed ReEDS runs
ReEDS_solve_year = 2035 #
Weather_year = 2012 #must be 2007-2013
Timestamps = 8760

pras_system = ReEDS2PRAS.reeds_to_pras(runpath,ReEDS_solve_year,NP,Timestamps,Weather_year) #returns a parameterized PRAS system
```

This will save out a pras system to the variable `pras_system` from the ReEDS2PRAS run.

## Acknowledgements
The developers are Luke Lavin and Surya Dhulipala. They acknowlege Gord Stephen for many helpfful comments for improving this package.