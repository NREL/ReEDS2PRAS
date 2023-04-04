using ReEDS2PRAS
using PRAS
using CSV
using DataFrames
using Dates
using HDF5
using Statistics
using Test
using TimeZones

include("testutils.jl")

@testset "ReEDS2PRAS" begin
    include("ntptests/ntpscenarios.jl")
end
