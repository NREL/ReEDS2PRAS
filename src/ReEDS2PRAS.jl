"""
Copyright 2022 Alliance for Sustainable Energy and other
NTP Project Developers. See the top-level COPYRIGHT file for details.
Author: Luke Lavin, Surya Chandan Dhulipala, Brandon Norton Benton, Hari Sundar
Email: luke.lavin@nrel.gov, suryachandan.dhulipala@nrel.gov, 
       brandon.benton@nrel.gov, hari.sundar@nrel.gov
"""
module ReEDS2PRAS

# Exports
export reeds_to_pras

# Imports
import CSV
import DataFrames
import Dates
import HDF5
import PRAS
import Statistics
import TimeZones
import InlineStrings
import JSON

# Includes
# Models
include("models/Region.jl")
include("models/Storage.jl")
include("models/Battery.jl")
include("models/Gen_Storage.jl")
include("models/Generator.jl")
include("models/Thermal_Gen.jl")
include("models/Variable_Gen.jl")
include("models/Line.jl")
include("models/utils.jl")

# Utils
include("utils/reeds_input_parsing.jl")
include("utils/runchecks.jl")
include("utils/reeds_data_parsing.jl")
include("utils/user_descriptors_parsing.jl")
#Main
include("main/extract_system_info.jl")
include("main/create_pras_system.jl")

# Module
include("reeds_to_pras.jl")

end
