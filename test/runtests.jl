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
    include("ntptests/ntpscenarios.jl")
end