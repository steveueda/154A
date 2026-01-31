function fuel = compute_range_gas(Range_req_km, V_outbound_mps, Preq_outbound_W, V_return_mps, Preq_return_W, assump)
% Computes fuel required for two-phase mission using BSFC (gasoline engine)
%
% Inputs:
%   Range_req_km - Required round trip range (km)
%   V_outbound_mps - Outbound cruise velocity (m/s)
%   Preq_outbound_W - Power required at outbound speed (W)
%   V_return_mps - Return cruise velocity (m/s) - best range speed
%   Preq_return_W - Power required at return speed (W)
%   assump - Assumptions struct with BSFC, fuel density, and margin parameters
%
% Outputs:
%   fuel.fuel_mass_req_kg - Required fuel mass (kg)
%   fuel.fuel_mass_total_kg - Total fuel mass with reserve (kg)
%   fuel.fuel_vol_total_L - Fuel volume (L)
%   fuel.tank_size_L - Tank size with margin (L)

Range_one_leg_km = Range_req_km / 2;
meters_per_km = 1000;
seconds_per_hour = 3600;
Range_one_leg_m = Range_one_leg_km * meters_per_km;

t_outbound_hr = Range_one_leg_m / (V_outbound_mps * seconds_per_hour);
t_return_hr = Range_one_leg_m / (V_return_mps * seconds_per_hour);

P_shaft_outbound_kW = (Preq_outbound_W / assump.eta_prop) / 1000;
P_shaft_return_kW = (Preq_return_W / assump.eta_prop) / 1000;

grams_per_kg = 1000;
BSFC_kg_per_kWh = assump.BSFC_g_per_kWh / grams_per_kg;
mdot_fuel_outbound_kg_per_hr = BSFC_kg_per_kWh * P_shaft_outbound_kW;
mdot_fuel_return_kg_per_hr = BSFC_kg_per_kWh * P_shaft_return_kW;

fuel_outbound_kg = mdot_fuel_outbound_kg_per_hr * t_outbound_hr;
fuel_return_kg = mdot_fuel_return_kg_per_hr * t_return_hr;
fuel_mass_required_kg = fuel_outbound_kg + fuel_return_kg;

fuel_mass_with_reserve_kg = fuel_mass_required_kg * (1 + assump.reserve_frac);
fuel_volume_L = fuel_mass_with_reserve_kg / assump.fuel_density_kg_per_L;
fuel_tank_size_L = fuel_volume_L * (1 + assump.tank_margin_frac);

fuel.fuel_mass_req_kg = fuel_mass_required_kg;
fuel.fuel_mass_total_kg = fuel_mass_with_reserve_kg;
fuel.fuel_vol_total_L = fuel_volume_L;
fuel.tank_size_L = fuel_tank_size_L;

end
