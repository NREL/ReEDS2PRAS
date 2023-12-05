using PRAS
using CSV
using DataFrames
using Dates
using HDF5
using Statistics
using Test
using TimeZones

using BenchmarkTools

using ReEDS2PRAS

# Use the main branch of R2P while running ReEDS 
# with the cases_repbenchmark.csv file. 
# update R2P path in the cases file before running ReEDS
#
# Then, run this file first with the main branch and 
# then the feature branch, record the reported mean time taken
# for both ReEDS2PRAS and PRAS
# while submitting pull request with the PR template

# Run ReEDS2PRAS

# Update with path to your ReEDS result
inpath = "path/to/reedsoutput"

ty = 2030

tp = joinpath(inpath, "v20231130_stress_WECC") 
# Hourly 2007-2013 = 7*8760 = 61320
WEATHERYEAR = 2007
n_timesteps = 61320

benchmark_r2p = @benchmark pras_sys = ReEDS2PRAS.reeds_to_pras(tp, ty, n_timesteps, WEATHERYEAR)



# Run PRAS

pras_sys = ReEDS2PRAS.reeds_to_pras(tp, ty, n_timesteps, WEATHERYEAR)

simulation = SequentialMonteCarlo(samples=10, seed=100)

benchmark_pras = @benchmark shortfall = assess(pras_sys,simulation,Shortfall())

println()
println("ReEDS2PRAS benchmark")

io = IOBuffer()
show(io, "text/plain", benchmark_r2p)
s = String(take!(io))
println(s)

println("\nPRAS benchmark with 10 samples")

io = IOBuffer()
show(io, "text/plain", benchmark_pras)
s = String(take!(io))
println(s)