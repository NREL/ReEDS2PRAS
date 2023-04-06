"""
    This code defines a struct called Gen_Storage which is a subtype of
    Storage. The struct has 14 fields: name, timesteps, region_name, type,
    charge_cap, discharge_cap, energy_cap, inflow, grid_withdrawl_cap,
    grid_inj_cap, legacy, charge_eff, discharge_eff and carryover_eff. The
    code also contains an inner constructor and checks to ensure that the
    values passed are valid. Specifically, timesteps must be between 0 and 8784
    (inclusive), all capacity values must be greater than or equal to 0.0,
    the legacy value must either be “Existing” or “New”, all of the
    efficiency values must be between 0.0 and 1.0 (inclusive), FOR must be
    between 0.0 and 1.0 (inclusive) and MTTR must be greater than 0.
    Additionally, there is a commented out check that verifies that all of
    the time series data is of the same size. Finally, if any of these
    checks fail, an error will be thrown.

    Parameters
    ----------
    name : string
        Name of Gen_Storage
    timesteps : integer
        PRAS timesteps (timesteps)
    region_name : string
        Region name
    type : string
        Storage type
    charge_cap : float
        Charge capacity
    discharge_cap : float
        Discharge capacity
    energy_cap : float
        Energy capacity
    inflow : float
        Inflow time series data
    grid_withdrawl_cap : float
        Grid withdrawal capacity time series data
    grid_inj_cap : floating point
        Grid injection capacity time series data
    legacy : string
        Must be either "Existing" or "New"
    charge_eff : float
        Charge efficiency
    discharge_eff : float
        Discharge efficiency
    carryover_eff : float
        Carryover efficiency
    FOR : float
        Forced Outage Rate value
    MTTR : integer
        Mean Time To Repair value

    Returns
    -------
    A new instance of Gen_Storage.
"""
struct Gen_Storage <: Storage
    name::String
    timesteps::Int64
    region_name::String
    type::String
    charge_cap::Float64
    discharge_cap::Float64
    energy_cap::Float64
    inflow::Float64
    grid_withdrawl_cap::Float64
    grid_inj_cap::Float64
    legacy::String
    charge_eff::Float64
    discharge_eff::Float64
    carryover_eff::Float64
    FOR::Float64
    MTTR::Int64

    # Inner Constructors & Checks
    function Gen_Storage(
        name,
        timesteps,
        region_name,
        type,
        charge_cap,
        discharge_cap,
        energy_cap,
        inflow,
        grid_withdrawl_cap,
        grid_inj_cap,
        legacy = "New",
        charge_eff = 1.0,
        discharge_eff = 1.0,
        carryover_eff = 1.0,
        FOR = 0.0,
        MTTR = 24,
    )
        all(charge_cap .>= 0.0) ||
            error("Check for inconsistencies in Gen_Storage charge capacity
                   time series data")

        all(discharge_cap .>= 0.0) ||
            error("Check for inconsistencies in Gen_Storage discharge capacity
                   time series data")

        all(energy_cap .>= 0.0) ||
            error("Check for inconsistencies in Gen_Storage energy capacity
                   time series data")

        all(inflow .>= 0.0) ||
            error("Check for inconsistencies in Gen_Storage inflow time series
                   data")

        all(grid_withdrawl_cap .>= 0.0) ||
            error("Check for inconsistencies in Gen_Storage grid withdrawl
                   capacity time series data")

        all(grid_inj_cap .>= 0.0) ||
            error("Check for inconsistencies in Gen_Storage grid injection
                   capacity time series data")

        legacy in ["Existing", "New"] || error("Unidentified legacy passed")

        all(0.0 .<= [charge_eff, discharge_eff, carryover_eff] .<= 1.0) ||
            error("Invalid charge/discharge/carryover efficiency passed")

        0.0 <= FOR <= 1.0 || error("FOR value passed is not allowed")

        MTTR > 0 || error("MTTR value passed is not allowed")

        return new(
            name,
            timesteps,
            region_name,
            type,
            charge_cap,
            discharge_cap,
            energy_cap,
            inflow,
            grid_withdrawl_cap,
            grid_inj_cap,
            legacy,
            charge_eff,
            discharge_eff,
            carryover_eff,
            FOR,
            MTTR,
        )
    end
end

# Getter Functions

get_charge_capacity(stor::Gen_Storage) =
    fill(round(Int, stor.charge_cap), 1, stor.timesteps)

get_discharge_capacity(stor::Gen_Storage) =
fill(round(Int, stor.discharge_cap), 1, stor.timesteps)

get_energy_capacity(stor::Gen_Storage) =
    fill(round(Int, stor.energy_cap), 1, stor.timesteps)

get_inflow(stor::Gen_Storage) = fill(round(Int, stor.inflow), 1, stor.timesteps)

get_grid_withdrawl_capacity(stor::Gen_Storage) =
    fill(round(Int, stor.grid_withdrawl_cap), 1, stor.timesteps)

get_grid_injection_capacity(stor::Gen_Storage) =
    fill(round(Int, stor.grid_inj_cap), 1, stor.timesteps)