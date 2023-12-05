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

# Set up this test
reedscase = "../reeds_cases/USA_VSC_2035"
solve_year = 2035
timesteps = 61320
weather_year = 2007
samples = 10
seed = 1

# Run ReEDS2PRAS
benchmark_r2p = @benchmark pras_sys = ReEDS2PRAS.reeds_to_pras(
    reedscase, solve_year, timesteps, weather_year)

# Run PRAS
pras_sys = ReEDS2PRAS.reeds_to_pras(
    reedscase, solve_year, timesteps, weather_year)

simulation = SequentialMonteCarlo(samples=samples, seed=seed)

benchmark_pras = @benchmark shortfall = assess(pras_sys, simulation, Shortfall())

println()
println("ReEDS2PRAS benchmark")

io = IOBuffer()
show(io, "text/plain", benchmark_r2p)
s = String(take!(io))
println(s)

println("\nPRAS benchmark with $(samples) samples")

io = IOBuffer()
show(io, "text/plain", benchmark_pras)
s = String(take!(io))
println(s)
