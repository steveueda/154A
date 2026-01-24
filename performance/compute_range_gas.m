function fuel = compute_range_gas(Range_req_km, V_cruise, Preq_W, assump)
% Computes fuel required for range using BSFC

t_req_hr = (Range_req_km * 1000) / V_cruise / 3600;
P_shaft_kW = (Preq_W / assump.eta_prop) / 1000;
mdot_fuel_kg_per_hr = (assump.BSFC_g_per_kWh / 1000) * P_shaft_kW;

fuel.fuel_mass_req_kg = mdot_fuel_kg_per_hr * t_req_hr;
fuel.fuel_mass_total_kg = fuel.fuel_mass_req_kg * (1 + assump.reserve_frac);
fuel.fuel_vol_total_L = fuel.fuel_mass_total_kg / assump.fuel_density_kg_per_L;
fuel.tank_size_L = fuel.fuel_vol_total_L * (1 + assump.tank_margin_frac);

end
