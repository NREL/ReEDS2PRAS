println("Running benchmark")

# Then, run this file first with the main branch and 
# then the feature branch, record the reported mean time taken 
# for both ReEDS2PRAS and PRAS, and the CONUS LOLE, nEUE output here 
# while submitting pull request with the PR template

# If making major changes to R2P model increase number of MC samples used

# Set up this test
reedscase = joinpath(dirname(rootfile),"reeds_cases","USA_VSC_2035")
solve_year = 2035
timesteps = 61320
weather_year = 2007
samples = 3
seed = 1

# Run ReEDS2PRAS
benchmark_r2p = @benchmark pras_sys = ReEDS2PRAS.reeds_to_pras(
    reedscase, solve_year, timesteps, weather_year)

# Run PRAS
pras_sys = ReEDS2PRAS.reeds_to_pras(
    reedscase, solve_year, timesteps, weather_year)

simulation = SequentialMonteCarlo(samples=samples, seed=seed)

benchmark_pras = @benchmark shortfall = assess(pras_sys, simulation, Shortfall())

shortfall = assess(pras_sys, simulation, Shortfall())

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

LOLE = PRAS.LOLE(shortfall[1]).lole.estimate
EUE = PRAS.EUE(shortfall[1]).eue.estimate
nEUE = EUE/sum(pras_sys.regions.load)

println("\nPRAS results")
println("LOLE=$(LOLE), nEUE=$(nEUE)")
