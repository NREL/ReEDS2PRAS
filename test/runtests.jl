using ReEDS2PRAS
using PRAS

using Dates
using CSV
using DataFrames
using Test
using HDF5

@testset "ReEDS2PRAS" begin
    # include("ercot/testlinesonly.jl")
    include("ercot/toyercot.jl")
end