function [W_final_N, S_final_m2, fuel_final_kg, tank_final_L, WeightHistory, converged, iter_converged, rel_error_final] = run_weight_iteration(req, assump, geom, fixedMasses, W0_N, components, env, dragParams)
% Weight convergence loop: weight -> wing sizing -> drag -> fuel -> new weight

aero = struct('ARw', geom.ARw, 'ew', assump.ew);
conv = unit_conversions();

W_old_N = W0_N;
iter_table = table();
converged = false;
iter_converged = assump.max_iter;
rel_error_final = 0;

for iter = 1:assump.max_iter
    Sref_m2 = size_wing(W_old_N, env.rho, req.V_stall_mps, assump.CL_max);
    
    components(1).Swet = Sref_m2 * 2;
    components(1).Lref = sqrt(Sref_m2 / geom.ARw);
    
    V_cruise = req.V_design_mps;
    drag = compute_drag(V_cruise, W_old_N, Sref_m2, components, env, aero, dragParams);
    power = compute_power(drag, V_cruise, assump);
    fuel = compute_range_gas(req.Range_req_km, V_cruise, power.Preq_W, assump);
    
    geomSI = struct('Sref_m2', Sref_m2, 'ARw', geom.ARw, 'taper_w', geom.taper_w, ...
        'tc_w', geom.tc_w, 'sweep_w_deg', geom.sweep_w_deg, 'lf_m', geom.lf_m, ...
        'wf_max_m', geom.wf_max_m, 'df_max_m', geom.df_max_m, 'Sh_m2', geom.Sh_m2, ...
        'ARh', geom.ARh, 'tc_h', geom.tc_h, 'Sv_m2', geom.Sv_m2, 'ARv', geom.ARv, ...
        'tc_v', geom.tc_v, 'Vmax_mps', assump.Vmax, 'n_ult', assump.n_ult);
    
    nicIn = buildNiccolaiInputs(geomSI, W_old_N, fixedMasses.engine_mass_kg);
    nicOut = niccolai_weight(nicIn);
    
    W_struct_kg = nicOut.W_struct_lb / conv.lb_per_kg;
    W_fixed_kg = fixedMasses.payload_mass_kg + fixedMasses.engine_mass_kg + ...
        fixedMasses.avionics_mass_kg + fixedMasses.landing_gear_mass_kg + fixedMasses.fuel_system_mass_kg;
    W_new_N = (W_struct_kg + W_fixed_kg + fuel.fuel_mass_total_kg) * conv.N_per_kg;
    
    rel_error = abs(W_new_N - W_old_N) / W_old_N;
    new_row = table(iter, W_old_N, Sref_m2, V_cruise, drag.CDp, drag.CDi, drag.CDtotal, ...
        power.Preq_W, fuel.fuel_mass_total_kg, W_new_N, rel_error, ...
        'VariableNames', {'iter', 'W_old_N', 'Sref_m2', 'V_cruise_used_mps', ...
        'CDp_at_cruise', 'CDi_at_cruise', 'CDtotal_at_cruise', ...
        'Preq_at_cruise_W', 'fuel_mass_total_kg', 'W_new_N', 'rel_error'});
    iter_table = [iter_table; new_row];
    
    if rel_error < assump.tol_frac
        converged = true;
        iter_converged = iter;
        rel_error_final = rel_error;
        break;
    end
    W_old_N = W_new_N;
    rel_error_final = rel_error;
end

W_final_N = W_new_N;
S_final_m2 = Sref_m2;
fuel_final_kg = fuel.fuel_mass_total_kg;
tank_final_L = fuel.tank_size_L;
WeightHistory = iter_table;

end
