using ReEDS2PRAS
using PRAS

using Dates
using Test
using HDF5

@testset "ReEDS2PRAS" begin
    include("ercot/toyercot.jl")
    include("ercot/testlinesonly.jl")
end