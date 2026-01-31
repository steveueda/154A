function [W_final_N, S_final_m2, fuel_final_kg, tank_final_L, WeightHistory, converged, iter_converged, rel_error_final] = run_weight_iteration(req, assump, geom, fixedMasses, W0_N, components, env, dragParams, k_corr)

aero = struct('ARw', geom.ARw, 'ew', assump.ew);
conv = unit_conversions();

W_old_N = W0_N;
iter_table = table();
converged = false;
iter_converged = assump.max_iter;
rel_error_final = 0;

wing_wet_area_factor = 2;
stall_margin_factor = 1.1;
stall_min_speed_mps = 10;
return_speed_limit_factor = 0.95;

payload_weight_N = fixedMasses.payload_mass_kg * conv.N_per_kg;

for iter = 1:assump.max_iter
    Sref_m2 = size_wing(W_old_N, env.rho, req.V_stall_mps, assump.CL_max);
    
    components(1).Swet = Sref_m2 * wing_wet_area_factor;
    components(1).Lref = sqrt(Sref_m2 / geom.ARw);
    
    W_outbound_N = W_old_N;
    W_return_N = W_old_N - payload_weight_N;
    
    V_outbound = req.V_design_mps;
    drag_outbound = compute_drag(V_outbound, W_outbound_N, Sref_m2, components, env, aero, dragParams);
    power_outbound = compute_power(drag_outbound, V_outbound, assump);
    
    V_stall_est_return = sqrt(2 * W_return_N / (env.rho * Sref_m2 * assump.CL_max));
    V_min_search = max(V_stall_est_return * stall_margin_factor, stall_min_speed_mps);
    V_max_search = V_outbound * return_speed_limit_factor;
    [V_return, Preq_return] = find_best_range_speed(W_return_N, Sref_m2, components, env, aero, dragParams, assump, V_min_search, V_max_search);
    
    if strcmpi(assump.prop_type, 'electric')
        fuel = compute_range_electric(req.Range_req_km, V_outbound, power_outbound.Preq_W, V_return, Preq_return, assump);
    else
        fuel = compute_range_gas(req.Range_req_km, V_outbound, power_outbound.Preq_W, V_return, Preq_return, assump);
    end
    
geomSI = struct('Sref_m2', Sref_m2, 'ARw', geom.ARw, 'taper_w', geom.taper_w, ...
    'tc_w', geom.tc_w, 'sweep_w_deg', geom.sweep_w_deg, 'lf_m', geom.lf_m, ...
    'wf_max_m', geom.wf_max_m, 'df_max_m', geom.df_max_m, 'Sh_m2', geom.Sh_m2, ...
    'ARh', geom.ARh, 'tc_h', geom.tc_h, 'Sv_m2', geom.Sv_m2, 'ARv', geom.ARv, ...
    'tc_v', geom.tc_v, 'Vmax_mps', req.V_design_mps, 'n_ult', assump.n_ult);
    
    nicIn = buildNiccolaiInputs(geomSI, W_outbound_N, fixedMasses.engine_mass_kg);
    nicOut = niccolai_weight(nicIn);
    
    W_wing_kg = (nicOut.W_wing_lb / conv.lb_per_kg) * k_corr.wing;
    W_fuse_kg = (nicOut.W_fuse_lb / conv.lb_per_kg) * k_corr.fuselage;
    W_ht_kg = (nicOut.W_ht_lb / conv.lb_per_kg) * k_corr.htail;
    W_vt_kg = (nicOut.W_vt_lb / conv.lb_per_kg) * k_corr.vtail;
    W_controls_kg = (nicOut.W_controls_lb / conv.lb_per_kg) * k_corr.controls;
    W_prop_install_kg = (nicOut.W_prop_install_lb / conv.lb_per_kg) * k_corr.prop_install;
    
    W_struct_kg = W_wing_kg + W_fuse_kg + W_ht_kg + W_vt_kg;
    W_fixed_kg = fixedMasses.payload_mass_kg + fixedMasses.engine_mass_kg + ...
        fixedMasses.avionics_mass_kg + fixedMasses.landing_gear_mass_kg + ...
        fixedMasses.fuel_system_mass_kg;
    W_new_N = (W_struct_kg + W_fixed_kg + fuel.fuel_mass_total_kg) * conv.N_per_kg;
    
    rel_error = abs(W_new_N - W_old_N) / W_old_N;
    new_row = table(iter, W_old_N, Sref_m2, V_outbound, drag_outbound.CDp, drag_outbound.CDi, drag_outbound.CDtotal, ...
        power_outbound.Preq_W, fuel.fuel_mass_total_kg, W_new_N, rel_error, ...
        'VariableNames', {'iter', 'W_old_N', 'Sref_m2', 'V_outbound_used_mps', ...
        'CDp_at_outbound', 'CDi_at_outbound', 'CDtotal_at_outbound', ...
        'Preq_at_outbound_W', 'fuel_mass_total_kg', 'W_new_N', 'rel_error'});
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
