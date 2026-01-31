function energy = compute_range_electric(Range_req_km, V_outbound_mps, Preq_outbound_W, V_return_mps, Preq_return_W, assump)
% Computes battery energy required for two-phase mission using electric propulsion
%
% Inputs:
%   Range_req_km - Required round trip range (km)
%   V_outbound_mps - Outbound cruise velocity (m/s)
%   Preq_outbound_W - Power required at outbound speed (W)
%   V_return_mps - Return cruise velocity (m/s) - best range speed
%   Preq_return_W - Power required at return speed (W)
%   assump - Assumptions struct with efficiency, energy density, and margin parameters
%
% Outputs:
%   energy.fuel_mass_req_kg - Required battery mass (kg)
%   energy.fuel_mass_total_kg - Total battery mass with reserve (kg)
%   energy.fuel_vol_total_L - Battery volume (L)
%   energy.tank_size_L - Battery size with margin (L)

Range_one_leg_km = Range_req_km / 2;
meters_per_km = 1000;
seconds_per_hour = 3600;
Range_one_leg_m = Range_one_leg_km * meters_per_km;

t_outbound_hr = Range_one_leg_m / (V_outbound_mps * seconds_per_hour);
t_return_hr = Range_one_leg_m / (V_return_mps * seconds_per_hour);

eta_total = assump.eta_prop * assump.eta_electric;
P_electric_outbound_W = Preq_outbound_W / eta_total;
P_electric_return_W = Preq_return_W / eta_total;

E_outbound_Wh = P_electric_outbound_W * t_outbound_hr;
E_return_Wh = P_electric_return_W * t_return_hr;
E_total_required_Wh = E_outbound_Wh + E_return_Wh;

battery_mass_required_kg = E_total_required_Wh / assump.battery_energy_density_Wh_per_kg;
battery_mass_with_reserve_kg = battery_mass_required_kg * (1 + assump.reserve_frac);

battery_density_kg_per_L = 1.7;
battery_volume_L = battery_mass_with_reserve_kg / battery_density_kg_per_L;
battery_size_with_margin_L = battery_volume_L * (1 + assump.tank_margin_frac);

energy.fuel_mass_req_kg = battery_mass_required_kg;
energy.fuel_mass_total_kg = battery_mass_with_reserve_kg;
energy.fuel_vol_total_L = battery_volume_L;
energy.tank_size_L = battery_size_with_margin_L;

end

