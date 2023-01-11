using ReEDS2PRAS
using PRAS
using Dates
using CSV
using DataFrames
using Test
using HDF5
using Statistics
using TimeZones

include("testutils.jl")

@testset "ReEDS2PRAS" begin
    # include("ercot/ntpscenarios.jl")
    # include("ercot/standardscenarios.jl")
    # include("ntp/ntp_scenarios.jl")
    include("ercot/toyercot.jl")
end