clear; clc;

script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'sizing'), fullfile(script_dir, 'weight'), ...
    fullfile(script_dir, 'aero'), fullfile(script_dir, 'performance'), ...
    fullfile(script_dir, 'stability'), fullfile(script_dir, 'utils'), ...
    fullfile(script_dir, 'report'));

CONV.mph_per_mps = 2.236936;
CONV.nmi_per_km = 1 / 1.852;
CONV.g_per_kg = 9.80665;
CONV.lb_per_kg = 2.20462262185;

run(fullfile(script_dir, 'user_inputs.m'));

assump = struct('CL_max', aero.CL_max, 'ew', aero.ew, ...
    'eta_prop', prop.eta_prop, ...
    'P_shaft_max_W', prop.P_shaft_max_W, 'prop_type', prop.type, ...
    'BSFC_g_per_kWh', prop.BSFC_g_per_kWh, ...
    'fuel_density_kg_per_L', prop.fuel_density_kg_per_L, ...
    'battery_energy_density_Wh_per_kg', prop.battery_energy_density_Wh_per_kg, ...
    'eta_electric', prop.eta_electric, ...
    'reserve_frac', prop.reserve_frac, ...
    'tank_margin_frac', prop.tank_margin_frac, 'n_ult', mass.n_ult, 'Q', aero.Q, ...
    'CD_misc', aero.CD_misc, 'CD_LP', aero.CD_LP, 'max_iter', ctrl.max_iter, ...
    'tol_frac', ctrl.tol_frac, ...
    'V_min_fraction', ctrl.V_min_fraction, 'V_max_fraction', ctrl.V_max_fraction);

components = struct('name', {}, 'type', {}, 'Swet', {}, 'Lref', {}, ...
    'tc', {}, 'xmc', {}, 'sweep_m_rad', {}, 'fineness', {});

wing_wet_area_factor = 2;
components(end+1) = struct('name',"wing",'type',"wing", 'Swet', NaN, 'Lref', NaN, ...
    'tc', geom.tc_w, 'xmc', geom.wing_xmc, 'sweep_m_rad', deg2rad(geom.sweep_w_deg), 'fineness', NaN);
components(end+1) = struct('name',"htail",'type',"wing", ...
    'Swet', geom.Sh_m2 * wing_wet_area_factor, 'Lref', sqrt(geom.Sh_m2 / geom.ARh), ...
    'tc', geom.tc_h, 'xmc', geom.htail_xmc, 'sweep_m_rad', deg2rad(geom.htail_sweep_deg), 'fineness', NaN);
components(end+1) = struct('name',"fuselage",'type',"fuselage", ...
    'Swet', geom.lf_m * (geom.wf_max_m + geom.df_max_m) * wing_wet_area_factor, 'Lref', geom.lf_m, ...
    'tc', NaN, 'xmc', NaN, 'sweep_m_rad', NaN, ...
    'fineness', geom.lf_m / max(geom.wf_max_m, geom.df_max_m));

env = struct('rho', geom.rho, 'mu', geom.mu, 'a', geom.a);
dragParams = struct('Q', aero.Q, 'CD_misc', aero.CD_misc, 'CD_LP', aero.CD_LP);
fixedMasses = struct('payload_mass_kg', req.payload_mass_kg, ...
    'engine_mass_kg', mass.engine_mass_kg, 'avionics_mass_kg', mass.avionics_mass_kg, ...
    'landing_gear_mass_kg', mass.landing_gear_mass_kg, 'fuel_system_mass_kg', mass.fuel_system_mass_kg);

[W_final_N, S_final_m2, fuel_final_kg, tank_final_L, WeightHistory, converged, iter_converged, rel_error_final] = ...
    run_weight_iteration(req, assump, geom, fixedMasses, mass.W0_N, components, env, dragParams, k_corr);

components(1).Swet = S_final_m2 * wing_wet_area_factor;
components(1).Lref = sqrt(S_final_m2 / geom.ARw);

payload_weight_N = fixedMasses.payload_mass_kg * CONV.g_per_kg;
W_return_N = W_final_N - payload_weight_N;

config = struct('W', W_final_N, 'W_return', W_return_N, 'Sref', S_final_m2, 'St', geom.Sh_m2, ...
    'ARw', geom.ARw, 'ARt', geom.ARh, 'components', components, ...
    'V_stall', req.V_stall_mps, 'V_cruise', req.V_design_mps, ...
    'fuel_mass_kg', fuel_final_kg, 'P_shaft_max', assump.P_shaft_max_W, ...
    'payload_mass_kg', fixedMasses.payload_mass_kg);

sweep_results_outbound = velocity_sweep(config, assump, env);

config_return = struct('W', W_return_N, 'W_return', W_return_N, 'Sref', S_final_m2, 'St', geom.Sh_m2, ...
    'ARw', geom.ARw, 'ARt', geom.ARh, 'components', components, ...
    'V_stall', req.V_stall_mps, 'V_cruise', req.V_design_mps, ...
    'fuel_mass_kg', fuel_final_kg, 'P_shaft_max', assump.P_shaft_max_W, ...
    'payload_mass_kg', 0);
sweep_results_return = velocity_sweep(config_return, assump, env);

sweep_results = sweep_results_outbound;

mission_range_backup = NaN;
try
    UserData = sweep_results.Properties.UserData;
    if ~isempty(UserData) && isstruct(UserData) && isfield(UserData, 'mission_range_nmi')
        mission_range_backup = UserData.mission_range_nmi;
        if isnan(mission_range_backup)
            warning('Mission range found in metadata but is NaN - check velocity_sweep calculation');
        end
    end
catch ME
    warning('Error accessing UserData: %s', ME.message);
end

requirements_metadata = struct('Range_req_nmi_round_trip', req.Range_req_km * CONV.nmi_per_km, ...
    'Top_speed_req_mph', req.V_design_mps * CONV.mph_per_mps, ...
    'Stall_speed_req_mph', req.V_stall_mps * CONV.mph_per_mps, ...
    'Payload_weight_lb', req.payload_mass_kg * CONV.lb_per_kg);

if isfield(sweep_results.Properties, 'UserData')
    sweep_results.Properties.UserData.requirements = requirements_metadata;
    if ~isnan(mission_range_backup)
        sweep_results.Properties.UserData.mission_range_nmi = mission_range_backup;
    end
else
    sweep_results.Properties.UserData = struct('requirements', requirements_metadata);
    if ~isnan(mission_range_backup)
        sweep_results.Properties.UserData.mission_range_nmi = mission_range_backup;
    end
end

V_vec_return = sweep_results_return.V_mps;
P_req_return = sweep_results_return.Preq_W;
P_avail_const_return = sweep_results_return.Pavail_W(1);

[~, idx_minP_return] = min(P_req_return);
[~, idx_maxROC_return] = max(sweep_results_return.ROC_mps);
[~, idx_maxRange_return] = max(sweep_results_return.Range_km);

V_intersect_list_return = [];
for i = 1:length(V_vec_return)-1
    if (P_req_return(i) - P_avail_const_return) * (P_req_return(i+1) - P_avail_const_return) <= 0
        V_int = V_vec_return(i) + (V_vec_return(i+1) - V_vec_return(i)) * ...
            (P_avail_const_return - P_req_return(i)) / (P_req_return(i+1) - P_req_return(i));
        V_intersect_list_return(end+1) = V_int;
    end
end

if ~isempty(V_intersect_list_return)
    V_max_return_mps = max(V_intersect_list_return);
    V_max_return_mph = V_max_return_mps * CONV.mph_per_mps;
else
    V_max_return_mps = NaN;
    V_max_return_mph = NaN;
end

CL_total_return = sweep_results_return.CL;
idx_stall_return = find(CL_total_return >= assump.CL_max, 1, 'last');
if isempty(idx_stall_return)
    [~, idx_stall_return] = min(abs(CL_total_return - assump.CL_max));
    V_stall_return_mps = V_vec_return(idx_stall_return);
elseif idx_stall_return == length(CL_total_return) || idx_stall_return == 1
    V_stall_return_mps = V_vec_return(idx_stall_return);
else
    V_stall_return_mps = V_vec_return(idx_stall_return) + (V_vec_return(idx_stall_return+1) - V_vec_return(idx_stall_return)) * ...
        (assump.CL_max - CL_total_return(idx_stall_return)) / (CL_total_return(idx_stall_return+1) - CL_total_return(idx_stall_return));
end
V_stall_return_mph = V_stall_return_mps * CONV.mph_per_mps;

design_speed_mph = 80;
V_80mph_mps = design_speed_mph / CONV.mph_per_mps;
if V_80mph_mps >= min(V_vec_return) && V_80mph_mps <= max(V_vec_return)
    Preq_80mph_return = interp1(V_vec_return, sweep_results_return.Preq_W, V_80mph_mps);
    ExcessP_80mph_return = interp1(V_vec_return, sweep_results_return.ExcessP_W, V_80mph_mps);
else
    Preq_80mph_return = NaN;
    ExcessP_80mph_return = NaN;
end

CL_sweep_return = sweep_results_return.CL;
CDtotal_sweep_return = sweep_results_return.CDtotal;
LD_sweep_return = CL_sweep_return ./ CDtotal_sweep_return;
[max_LD_return, idx_maxLD_return] = max(LD_sweep_return);
V_maxLD_return_mps = V_vec_return(idx_maxLD_return);
V_maxLD_return_mph = V_maxLD_return_mps * CONV.mph_per_mps;


lambda = geom.taper_w;
b = sqrt(geom.ARw * S_final_m2);
c_root = (2 * S_final_m2) / (b * (1 + lambda));
c_tip = lambda * c_root;
mac_coeff = 2/3;
cbar = mac_coeff * c_root * (1 + lambda + lambda^2) / (1 + lambda);
aerodynamic_center_fraction = 0.25;

stabilityInputs = struct('cbar_m', cbar, 'x_ac_w_m', aerodynamic_center_fraction*cbar, ...
    'Sh_m2', geom.Sh_m2, 'eta_tail', stab.eta_tail, 'a_w', stab.a_w, 'a_t', stab.a_t, ...
    'deps_dalpha', stab.deps_dalpha, 'lt_m', stab.lt_frac*geom.lf_m);

compMassLoc_outbound = struct('name', {}, 'mass_kg', {}, 'x_m', {});
compMassLoc_outbound(end+1) = struct('name', "payload", 'mass_kg', fixedMasses.payload_mass_kg, 'x_m', cg.payload_x_m);
compMassLoc_outbound(end+1) = struct('name', "engine", 'mass_kg', fixedMasses.engine_mass_kg, 'x_m', cg.engine_x_m);
compMassLoc_outbound(end+1) = struct('name', "avionics", 'mass_kg', fixedMasses.avionics_mass_kg, 'x_m', cg.avionics_x_m);
compMassLoc_outbound(end+1) = struct('name', "fuel_system", 'mass_kg', fixedMasses.fuel_system_mass_kg, 'x_m', cg.fuel_system_x_m);
compMassLoc_outbound(end+1) = struct('name', "landing_gear", 'mass_kg', fixedMasses.landing_gear_mass_kg, 'x_m', cg.landing_gear_x_m);

if strcmpi(prop.type, 'electric')
    compMassLoc_outbound(end+1) = struct('name', "battery", 'mass_kg', fuel_final_kg, 'x_m', cg.battery_x_m);
else
    compMassLoc_outbound(end+1) = struct('name', "fuel", 'mass_kg', fuel_final_kg, 'x_m', cg.fuel_x_m);
end

compMassLoc_return = struct('name', {}, 'mass_kg', {}, 'x_m', {});
compMassLoc_return(end+1) = struct('name', "engine", 'mass_kg', fixedMasses.engine_mass_kg, 'x_m', cg.engine_x_m);
compMassLoc_return(end+1) = struct('name', "avionics", 'mass_kg', fixedMasses.avionics_mass_kg, 'x_m', cg.avionics_x_m);
compMassLoc_return(end+1) = struct('name', "fuel_system", 'mass_kg', fixedMasses.fuel_system_mass_kg, 'x_m', cg.fuel_system_x_m);
compMassLoc_return(end+1) = struct('name', "landing_gear", 'mass_kg', fixedMasses.landing_gear_mass_kg, 'x_m', cg.landing_gear_x_m);

if strcmpi(prop.type, 'electric')
    compMassLoc_return(end+1) = struct('name', "battery", 'mass_kg', fuel_final_kg, 'x_m', cg.battery_x_m);
else
    compMassLoc_return(end+1) = struct('name', "fuel", 'mass_kg', fuel_final_kg, 'x_m', cg.fuel_x_m);
end

[x_cg_outbound_m, total_mass_outbound_kg, CGTable_outbound] = compute_cg(compMassLoc_outbound);
[x_cg_return_m, total_mass_return_kg, CGTable_return] = compute_cg(compMassLoc_return);
x_np_from_wing_LE = compute_neutral_point(stabilityInputs, S_final_m2, geom.ARw);

x_cg_outbound_from_wing_LE = x_cg_outbound_m - cg.wing_LE_from_front_m;
x_cg_return_from_wing_LE = x_cg_return_m - cg.wing_LE_from_front_m;
[SM_outbound, ~] = compute_static_margin(x_cg_outbound_from_wing_LE, x_np_from_wing_LE, cbar, stab.SM_min, stab.SM_max);
[SM_return, ~] = compute_static_margin(x_cg_return_from_wing_LE, x_np_from_wing_LE, cbar, stab.SM_min, stab.SM_max);
x_np_m = x_np_from_wing_LE + cg.wing_LE_from_front_m;

x_cg_m = x_cg_outbound_m;
CGTable = CGTable_outbound;
SM = SM_outbound;

geomSI_final = struct('Sref_m2', S_final_m2, 'ARw', geom.ARw, 'taper_w', geom.taper_w, ...
    'tc_w', geom.tc_w, 'sweep_w_deg', geom.sweep_w_deg, 'lf_m', geom.lf_m, ...
    'wf_max_m', geom.wf_max_m, 'df_max_m', geom.df_max_m, 'Sh_m2', geom.Sh_m2, ...
    'ARh', geom.ARh, 'tc_h', geom.tc_h, 'Sv_m2', geom.Sv_m2, 'ARv', geom.ARv, ...
    'tc_v', geom.tc_v, 'Vmax_mps', req.V_design_mps, 'n_ult', assump.n_ult);
nicIn_final = buildNiccolaiInputs(geomSI_final, W_final_N, fixedMasses.engine_mass_kg);
nicOut_final = niccolai_weight(nicIn_final);
conv = unit_conversions();

V_vec = sweep_results.V_mps;
P_req = sweep_results.Preq_W;
P_avail_const = sweep_results.Pavail_W(1);

[~, idx_minP] = min(P_req);
[~, idx_maxROC] = max(sweep_results.ROC_mps);
[~, idx_maxRange] = max(sweep_results.Range_km);

V_intersect_list = [];
for i = 1:length(V_vec)-1
    if (P_req(i) - P_avail_const) * (P_req(i+1) - P_avail_const) <= 0
        V_int = V_vec(i) + (V_vec(i+1) - V_vec(i)) * ...
            (P_avail_const - P_req(i)) / (P_req(i+1) - P_req(i));
        V_intersect_list(end+1) = V_int;
    end
end

if ~isempty(V_intersect_list)
    V_max_mps = max(V_intersect_list);
    V_max_mph = V_max_mps * CONV.mph_per_mps;
else
    V_max_mps = NaN;
    V_max_mph = NaN;
end

CL_total = sweep_results.CL;
idx_stall = find(CL_total >= assump.CL_max, 1, 'last');
if isempty(idx_stall)
    [~, idx_stall] = min(abs(CL_total - assump.CL_max));
    V_stall_mps = V_vec(idx_stall);
elseif idx_stall == length(CL_total) || idx_stall == 1
    V_stall_mps = V_vec(idx_stall);
else
    V_stall_mps = V_vec(idx_stall) + (V_vec(idx_stall+1) - V_vec(idx_stall)) * ...
        (assump.CL_max - CL_total(idx_stall)) / (CL_total(idx_stall+1) - CL_total(idx_stall));
end
V_stall_mph = V_stall_mps * CONV.mph_per_mps;

design_speed_mph = 80;
V_80mph_mps = design_speed_mph / CONV.mph_per_mps;
if V_80mph_mps >= min(V_vec) && V_80mph_mps <= max(V_vec)
    Preq_80mph = interp1(V_vec, sweep_results.Preq_W, V_80mph_mps);
    ExcessP_80mph = interp1(V_vec, sweep_results.ExcessP_W, V_80mph_mps);
else
    V_80mph_mps = NaN;
    Preq_80mph = NaN;
    ExcessP_80mph = NaN;
end
results = struct();
results.converged = converged;
results.iter_converged = iter_converged;
results.rel_error_final = rel_error_final;
results.max_iter = assump.max_iter;
results.W_final_N = W_final_N;
results.W_final_kg = W_final_N / CONV.g_per_kg;
results.S_final_m2 = S_final_m2;
results.WeightHistory = WeightHistory;
results.Range_req_km = req.Range_req_km;
results.V_top_req_mph = req.V_design_mps * CONV.mph_per_mps;
results.V_stall_req_mph = req.V_stall_mps * CONV.mph_per_mps;

try
    UserData = sweep_results.Properties.UserData;
    if ~isempty(UserData) && isstruct(UserData) && isfield(UserData, 'mission_range_nmi')
        results.Range_mission_nmi = UserData.mission_range_nmi;
        results.Range_mission_km = results.Range_mission_nmi / CONV.nmi_per_km;
        if isnan(results.Range_mission_nmi)
            warning('Mission range is NaN - check velocity_sweep calculation');
        end
    else
        results.Range_mission_nmi = NaN;
        results.Range_mission_km = NaN;
        warning('Mission range not found in sweep_results metadata');
    end
catch ME
    results.Range_mission_nmi = NaN;
    results.Range_mission_km = NaN;
    warning('Mission range not found in sweep_results metadata: %s', ME.message);
end
results.Range_best_km = sweep_results.Range_km(idx_maxRange);
results.Range_best_nmi = results.Range_best_km * CONV.nmi_per_km;
results.V_max_mph = V_max_mph;
results.V_max_mps = V_max_mps;
results.V_stall_mph = V_stall_mph;
results.V_stall_mps = V_stall_mps;
results.V_minPower_mps = V_vec(idx_minP);
results.V_minPower_mph = results.V_minPower_mps * CONV.mph_per_mps;
results.Preq_min_W = P_req(idx_minP);
results.V_bestRange_mps = V_vec(idx_maxRange);
results.V_bestRange_mph = results.V_bestRange_mps * CONV.mph_per_mps;
results.V_maxROC_mps = V_vec(idx_maxROC);
results.V_maxROC_mph = results.V_maxROC_mps * CONV.mph_per_mps;
results.ROC_max_mps = sweep_results.ROC_mps(idx_maxROC);
results.gamma_maxROC_deg = sweep_results.gamma_deg(idx_maxROC);
results.V_80mph_mps = V_80mph_mps;
results.Preq_80mph_W = Preq_80mph;
results.ExcessP_80mph_W = ExcessP_80mph;
percent_conversion = 100;
results.xCG_outbound_m = x_cg_outbound_m;
results.xCG_return_m = x_cg_return_m;
results.xCG_m = x_cg_outbound_m;
results.xNP_m = x_np_m;
results.SM_outbound = SM_outbound;
results.SM_return = SM_return;
results.SM = SM_outbound;
results.CG_outbound_pctMAC = percent_conversion * x_cg_outbound_from_wing_LE / cbar;
results.CG_return_pctMAC = percent_conversion * x_cg_return_from_wing_LE / cbar;
results.CG_pctMAC = percent_conversion * x_cg_outbound_from_wing_LE / cbar;
results.NP_pctMAC = percent_conversion * x_np_from_wing_LE / cbar;
results.SM_outbound_pctMAC = percent_conversion * SM_outbound;
results.SM_return_pctMAC = percent_conversion * SM_return;
results.SM_pctMAC = percent_conversion * SM_outbound;
results.CGTable_outbound = CGTable_outbound;
results.CGTable_return = CGTable_return;
results.CGTable = CGTable;
results.W_outbound_kg = total_mass_outbound_kg;
results.W_return_kg = total_mass_return_kg;
results.sweep_results = sweep_results;
results.sweep_results_outbound = sweep_results_outbound;
results.sweep_results_return = sweep_results_return;
results.V_max_return_mps = V_max_return_mps;
results.V_max_return_mph = V_max_return_mph;
results.V_stall_return_mps = V_stall_return_mps;
results.V_stall_return_mph = V_stall_return_mph;
results.V_minPower_return_mps = V_vec_return(idx_minP_return);
results.V_minPower_return_mph = results.V_minPower_return_mps * CONV.mph_per_mps;
results.Preq_min_return_W = P_req_return(idx_minP_return);
results.V_bestRange_return_mps = V_vec_return(idx_maxRange_return);
results.V_bestRange_return_mph = results.V_bestRange_return_mps * CONV.mph_per_mps;
results.Range_best_return_km = sweep_results_return.Range_km(idx_maxRange_return);
results.Range_best_return_nmi = results.Range_best_return_km * CONV.nmi_per_km;
results.V_maxROC_return_mps = V_vec_return(idx_maxROC_return);
results.V_maxROC_return_mph = results.V_maxROC_return_mps * CONV.mph_per_mps;
results.ROC_max_return_mps = sweep_results_return.ROC_mps(idx_maxROC_return);
results.gamma_maxROC_return_deg = sweep_results_return.gamma_deg(idx_maxROC_return);
results.Preq_80mph_return_W = Preq_80mph_return;
results.ExcessP_80mph_return_W = ExcessP_80mph_return;
results.max_LD_return = max_LD_return;
results.V_maxLD_return_mps = V_maxLD_return_mps;
results.V_maxLD_return_mph = V_maxLD_return_mph;
results.W_wing_kg = (nicOut_final.W_wing_lb / conv.lb_per_kg) * k_corr.wing;
results.W_fuse_kg = (nicOut_final.W_fuse_lb / conv.lb_per_kg) * k_corr.fuselage;
results.W_ht_kg = (nicOut_final.W_ht_lb / conv.lb_per_kg) * k_corr.htail;
results.W_vt_kg = (nicOut_final.W_vt_lb / conv.lb_per_kg) * k_corr.vtail;
results.W_struct_kg = results.W_wing_kg + results.W_fuse_kg + results.W_ht_kg + results.W_vt_kg;
results.W_controls_kg = (nicOut_final.W_controls_lb / conv.lb_per_kg) * k_corr.controls;
results.W_prop_install_kg = (nicOut_final.W_prop_install_lb / conv.lb_per_kg) * k_corr.prop_install;
results.fuel_final_kg = fuel_final_kg;

printResults(results);

D_parasite = sweep_results.Drag_N .* (sweep_results.CDp ./ sweep_results.CDtotal);
D_induced = sweep_results.Drag_N .* (sweep_results.CDi ./ sweep_results.CDtotal);
P_parasite = D_parasite .* V_vec;
P_induced = D_induced .* V_vec;
P_avail_vec = sweep_results.Pavail_W;

V_intersect_plot = [];
for i = 1:length(V_vec)-1
    if (P_req(i) - P_avail_vec(i)) * (P_req(i+1) - P_avail_vec(i+1)) <= 0
        V_int = V_vec(i) + (V_vec(i+1) - V_vec(i)) * ...
            (P_avail_vec(i) - P_req(i)) / (P_req(i+1) - P_req(i) - P_avail_vec(i+1) + P_avail_vec(i));
        V_intersect_plot(end+1) = V_int;
    end
end

Pmin = P_req(idx_minP);
V_minP = V_vec(idx_minP);
V_Range_max = V_vec(idx_maxRange);
P_rng = P_req(idx_maxRange);

ylim_margin = 1.1;
marker_size = 4;
font_size = 11;
title_font_size = 14;
label_font_size = 12;

% Power plot with dark mode styling
fig1 = figure('Color', 'k');
ax1 = axes('Parent', fig1, 'Color', 'k', 'XColor', [0.85 0.85 0.85], ...
    'YColor', [0.85 0.85 0.85], 'GridColor', [0.3 0.3 0.3], 'GridAlpha', 0.5);
hold(ax1, 'on');
grid(ax1, 'on');
grid(ax1, 'minor');

% Shade excess power region (where P_avail > P_req)
% Fill the area between P_req and P_avail curves
% The excess power region (where P_avail > P_req) will be clearly visible
fill(ax1, [V_vec; flipud(V_vec)], [P_req; flipud(P_avail_vec)], ...
    [0.1 0.3 0.6], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% Main curves with different shades of blue for dark mode
plot(ax1, V_vec, P_req, 'Color', [0.3 0.7 1.0], 'LineWidth', 2.5, 'DisplayName', 'Power Required');
plot(ax1, V_vec, P_avail_vec, 'Color', [0.1 0.4 0.8], 'LineWidth', 2.5, 'DisplayName', 'Power Available');

% Calculate offsets for text labels to avoid curve intersections (used when coordinates are NaN)
y_range = ylim_margin * max(P_avail_vec);
y_offset = 0.03 * y_range;  % 3% of y-axis range for vertical offset
x_offset = 0.02 * (max(V_vec) - min(V_vec));  % 2% of x-axis range for horizontal offset

% Label stall speed
P_stall = interp1(V_vec, P_req, V_stall_mps);
plot(ax1, [V_stall_mps V_stall_mps], [0 P_stall], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5, 'HandleVisibility', 'off');
plot(ax1, V_stall_mps, P_stall, 'o', 'MarkerSize', marker_size, 'MarkerFaceColor', [1 1 1], ...
    'MarkerEdgeColor', [1 1 1], 'LineWidth', 1.5, 'HandleVisibility', 'off');
if isnan(ctrl.label_power.v_stall_x)
    label_x_stall = V_stall_mps + 0.5*x_offset;
else
    label_x_stall = ctrl.label_power.v_stall_x;
end
if isnan(ctrl.label_power.v_stall_y)
    label_y_stall = P_stall + 0.5*y_offset;
else
    label_y_stall = ctrl.label_power.v_stall_y;
end
text(ax1, label_x_stall, label_y_stall, 'v_{stall}', ...
    'VerticalAlignment', 'bottom', 'FontSize', font_size, 'Color', [0.9 0.9 0.9]);

% Label minimum power speed
plot(ax1, [V_minP V_minP], [0 Pmin], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5, 'HandleVisibility', 'off');
plot(ax1, V_minP, Pmin, 'o', 'MarkerSize', marker_size, 'MarkerFaceColor', [1 1 1], ...
    'MarkerEdgeColor', [1 1 1], 'LineWidth', 1.5, 'HandleVisibility', 'off');
if isnan(ctrl.label_power.v_minpower_x)
    label_x_minpower = V_minP - 0*x_offset;
else
    label_x_minpower = ctrl.label_power.v_minpower_x;
end
if isnan(ctrl.label_power.v_minpower_y)
    label_y_minpower = Pmin - 2.5*y_offset;
else
    label_y_minpower = ctrl.label_power.v_minpower_y;
end
text(ax1, label_x_minpower, label_y_minpower, 'v_{minpower}', ...
    'VerticalAlignment', 'top', 'FontSize', font_size, 'Color', [0.9 0.9 0.9]);

% Label best range speed
plot(ax1, [V_Range_max V_Range_max], [0 P_rng], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5, 'HandleVisibility', 'off');
plot(ax1, V_Range_max, P_rng, 'o', 'MarkerSize', marker_size, 'MarkerFaceColor', [1 1 1], ...
    'MarkerEdgeColor', [1 1 1], 'LineWidth', 1.5, 'HandleVisibility', 'off');
if isnan(ctrl.label_power.v_range_x)
    label_x_range = V_Range_max - 0.2*x_offset;
else
    label_x_range = ctrl.label_power.v_range_x;
end
if isnan(ctrl.label_power.v_range_y)
    label_y_range = P_rng + y_offset;
else
    label_y_range = ctrl.label_power.v_range_y;
end
text(ax1, label_x_range, label_y_range, 'v_{range}', ...
    'VerticalAlignment', 'bottom', 'FontSize', font_size, 'Color', [0.9 0.9 0.9]);

% Label maximum speed
if ~isempty(V_intersect_plot) && ~isnan(V_max_mps)
    Vmax = V_max_mps;
    Pmax = interp1(V_vec, P_req, Vmax);
    plot(ax1, [Vmax Vmax], [0 Pmax], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    plot(ax1, Vmax, Pmax, 'o', 'MarkerSize', marker_size, 'MarkerFaceColor', [1 1 1], ...
        'MarkerEdgeColor', [1 1 1], 'LineWidth', 1.5, 'HandleVisibility', 'off');
    if isnan(ctrl.label_power.v_max_x)
        label_x_max = Vmax + 1.5*x_offset;
    else
        label_x_max = ctrl.label_power.v_max_x;
    end
    if isnan(ctrl.label_power.v_max_y)
        label_y_max = Pmax + 0.5*y_offset;
    else
        label_y_max = ctrl.label_power.v_max_y;
    end
    text(ax1, label_x_max, label_y_max, 'v_{max}', ...
        'VerticalAlignment', 'bottom', 'FontSize', font_size, 'Color', [0.9 0.9 0.9]);
end

xlabel(ax1, 'Velocity (m/s)', 'FontSize', label_font_size, 'Color', [0.9 0.9 0.9]);
ylabel(ax1, 'Power (W)', 'FontSize', label_font_size, 'Color', [0.9 0.9 0.9]);
title(ax1, 'Power vs. Velocity', 'FontSize', title_font_size, 'Color', [0.95 0.95 0.95]);
xlim(ax1, [ctrl.power_plot_xlim_low*min(V_vec) ctrl.power_plot_xlim_high*max(V_vec)]);
ylim(ax1, [0 ylim_margin*max(P_avail_vec)]);
legend(ax1, 'Location', 'best', 'TextColor', [0.9 0.9 0.9], 'Color', [0.15 0.15 0.15], 'EdgeColor', [0.4 0.4 0.4]);

% Drag plot with dark mode styling
fig2 = figure('Color', 'k');
ax2 = axes('Parent', fig2, 'Color', 'k', 'XColor', [0.85 0.85 0.85], ...
    'YColor', [0.85 0.85 0.85], 'GridColor', [0.3 0.3 0.3], 'GridAlpha', 0.5);
hold(ax2, 'on');
grid(ax2, 'on');
grid(ax2, 'minor');

Drag_total = sweep_results.Drag_N;
plot(ax2, V_vec, Drag_total, 'Color', [0.3 0.7 1.0], 'LineWidth', 2.5, 'DisplayName', 'Total Drag');
plot(ax2, V_vec, D_parasite, 'Color', [0.1 0.4 0.8], 'LineWidth', 1.5, 'LineStyle', '-', 'DisplayName', 'Parasite Drag');
plot(ax2, V_vec, D_induced, 'Color', [0.5 0.8 1.0], 'LineWidth', 1.5, 'LineStyle', '-', 'DisplayName', 'Induced Drag');

% Find minimum drag speed
[~, idx_minDrag] = min(Drag_total);
V_minDrag = V_vec(idx_minDrag);
Drag_minDrag = Drag_total(idx_minDrag);

% Calculate offsets for text labels to avoid curve intersections (used when coordinates are NaN)
y_range_drag = ylim_margin * max(Drag_total);
y_offset_drag = 0.03 * y_range_drag;  % 3% of y-axis range for vertical offset
x_offset_drag = 0.02 * (max(V_vec) - min(V_vec));  % 2% of x-axis range for horizontal offset

% Label minimum drag speed
plot(ax2, [V_minDrag V_minDrag], [0 Drag_minDrag], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5, 'HandleVisibility', 'off');
plot(ax2, V_minDrag, Drag_minDrag, 'o', 'MarkerSize', marker_size, 'MarkerFaceColor', [1 1 1], ...
    'MarkerEdgeColor', [1 1 1], 'LineWidth', 1.5, 'HandleVisibility', 'off');
if isnan(ctrl.label_drag.v_mindrag_x)
    label_x_mindrag = V_minDrag + 0.3*x_offset_drag;
else
    label_x_mindrag = ctrl.label_drag.v_mindrag_x;
end
if isnan(ctrl.label_drag.v_mindrag_y)
    label_y_mindrag = Drag_minDrag + y_offset_drag;
else
    label_y_mindrag = ctrl.label_drag.v_mindrag_y;
end
text(ax2, label_x_mindrag, label_y_mindrag, 'v_{mindrag}', ...
    'VerticalAlignment', 'bottom', 'FontSize', font_size, 'Color', [0.9 0.9 0.9]);

% Label stall speed
Drag_stall = interp1(V_vec, Drag_total, V_stall_mps);
plot(ax2, [V_stall_mps V_stall_mps], [0 Drag_stall], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5, 'HandleVisibility', 'off');
plot(ax2, V_stall_mps, Drag_stall, 'o', 'MarkerSize', marker_size, ...
    'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [1 1 1], 'LineWidth', 1.5, 'HandleVisibility', 'off');
if isnan(ctrl.label_drag.v_stall_x)
    label_x_stall_drag = V_stall_mps + 0.3*x_offset_drag;
else
    label_x_stall_drag = ctrl.label_drag.v_stall_x;
end
if isnan(ctrl.label_drag.v_stall_y)
    label_y_stall_drag = Drag_stall - 2.0*y_offset_drag;
else
    label_y_stall_drag = ctrl.label_drag.v_stall_y;
end
text(ax2, label_x_stall_drag, label_y_stall_drag, 'v_{stall}', ...
    'VerticalAlignment', 'top', 'FontSize', font_size, 'Color', [0.9 0.9 0.9]);

xlabel(ax2, 'Velocity (m/s)', 'FontSize', label_font_size, 'Color', [0.9 0.9 0.9]);
ylabel(ax2, 'Drag (N)', 'FontSize', label_font_size, 'Color', [0.9 0.9 0.9]);
title(ax2, 'Drag vs Velocity', 'FontSize', title_font_size, 'Color', [0.95 0.95 0.95]);
xlim(ax2, [ctrl.drag_plot_xlim_low*min(V_vec) ctrl.drag_plot_xlim_high*max(V_vec)]);
ylim(ax2, [0 ylim_margin*max(Drag_total)]);
legend(ax2, 'Location', 'best', 'TextColor', [0.9 0.9 0.9], 'Color', [0.15 0.15 0.15], 'EdgeColor', [0.4 0.4 0.4]);
