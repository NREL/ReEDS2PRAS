"""
    This code defines a struct called Battery which is a subtype of
    Storage. The struct has 13 fields: name, timesteps, region_name, type,
    charge_cap, discharge_cap, energy_cap, legacy, charge_eff,
    discharge_eff, carryover_eff, FOR and MTTR. The code also includes an
    inner constructor and checks to ensure that the values passed are
    valid. The constructor checks that timesteps is between 0 and 8784
    (inclusive), charge_cap and discharge_cap are greater than 0.0, energy_cap
    is greater than 0.0, legacy is either "Existing" or "New", all of the
    efficiency values are between 0.0 and 1.0 (inclusive), FOR is between 0.0
    and 1.0 (inclusive) and MTTR is greater than 0. If any of these checks fail
    an error will be thrown.

    Parameters
    ----------
    name : String
        The name of the battery
    timesteps : Int64
        Number of PRAS timesteps
    region_name: String
        Name of the region
    type : String
        Type of battery
    charge_cap : Float64
        Charge capacity
    discharge_cap : Float64
        Discharge capacity
    energy_cap : Float64
        Energy capacity
    legacy : String
        Battery's legacy (existing or new)
    charge_eff : Float64
        Charge efficiency
    discharge_eff : Float64
        Discharge efficiency
    carryover_eff : Float64
        Carryover efficiency
    FOR : Float64
        Factor of restoration
    MTTR : Int64
        Mean time to restore

    Returns
    -------
    Struct with properties related to the batter
"""
struct Battery <: Storage
    name::String
    timesteps::Int64
    region_name::String
    type::String
    charge_cap::Float64
    discharge_cap::Float64
    energy_cap::Float64
    legacy::String
    charge_eff::Float64
    discharge_eff::Float64
    carryover_eff::Float64
    FOR::Float64
    MTTR::Int64

    # Inner Constructors & Checks
    function Battery(
        name,
        timesteps,
        region_name,
        type,
        charge_cap,
        discharge_cap,
        energy_cap,
        legacy = "New",
        charge_eff = 1.0,
        discharge_eff = 1.0,
        carryover_eff = 1.0,
        FOR = 0.0,
        MTTR = 24,
    )
        @debug "cap_P = $(discharge_cap) MW and cap_E = $(energy_cap) MWh"

        charge_cap > 0.0 || error("Charge capacity passed is not allowed")

        discharge_cap > 0.0 || error("Discharge capacity passed is not allowed")

        energy_cap > 0.0 || error(
            "Energy capacity passed is not allowed: " *
            "$(name) $(charge_cap)MW $(energy_cap)MWh",
        )

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

get_charge_capacity(stor::Battery) = fill(round(Int, stor.charge_cap), 1, stor.timesteps)

get_discharge_capacity(stor::Battery) =
    fill(round(Int, stor.discharge_cap), 1, stor.timesteps)

get_energy_capacity(stor::Battery) = fill(round(Int, stor.energy_cap), 1, stor.timesteps)