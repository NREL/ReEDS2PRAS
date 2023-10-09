# Revert changes here before merging
#using Revise

using ReEDS2PRAS
using PRAS
using CSV
using DataFrames
using Dates
using HDF5
using Statistics
using Test
using TimeZones

#using Logging 

include("testutils.jl")

#io = open("log.txt","w+")
#logger = ConsoleLogger(io,Logging.Debug)

#global_logger(logger)

@testset "ReEDS2PRAS" begin
    include("ntptests/ntpscenarios.jl")
end

#close(io)
