%% MAIN REQUIREMENTS ITERATION

clear; clc;

%% SETUP
script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'sizing'), fullfile(script_dir, 'weight'), ...
    fullfile(script_dir, 'aero'), fullfile(script_dir, 'performance'), ...
    fullfile(script_dir, 'stability'), fullfile(script_dir, 'utils'), ...
    fullfile(script_dir, 'report'));

CONV.mph_per_mps = 2.236936;
CONV.nmi_per_km = 1 / 1.852;
CONV.g_per_kg = 9.80665;

%% LOAD INPUTS
run(fullfile(script_dir, 'user_inputs.m'));

assump = struct('CL_max', aero.CL_max, 'ew', aero.ew, ...
    'eta_prop', prop.eta_prop, ...
    'P_shaft_max_W', prop.P_shaft_max_W, 'BSFC_g_per_kWh', prop.BSFC_g_per_kWh, ...
    'fuel_density_kg_per_L', prop.fuel_density_kg_per_L, 'reserve_frac', prop.reserve_frac, ...
    'tank_margin_frac', prop.tank_margin_frac, 'n_ult', mass.n_ult, 'Q', aero.Q, ...
    'CD_misc', aero.CD_misc, 'CD_LP', aero.CD_LP, 'max_iter', ctrl.max_iter, ...
    'tol_frac', ctrl.tol_frac, 'Vmin', ctrl.Vmin, 'Vmax', ctrl.Vmax, 'N_sweep', ctrl.N_sweep);

%% COMPONENT DEFINITIONS FOR DRAG
components = struct('name', {}, 'type', {}, 'Swet', {}, 'Lref', {}, ...
    'tc', {}, 'xmc', {}, 'sweep_m_rad', {}, 'fineness', {});

components(end+1) = struct('name',"wing",'type',"wing", 'Swet', NaN, 'Lref', NaN, ...
    'tc', geom.tc_w, 'xmc', geom.wing_xmc, 'sweep_m_rad', deg2rad(geom.sweep_w_deg), 'fineness', NaN);
components(end+1) = struct('name',"htail",'type',"wing", ...
    'Swet', geom.Sh_m2 * 2, 'Lref', sqrt(geom.Sh_m2 / geom.ARh), ...
    'tc', geom.tc_h, 'xmc', geom.htail_xmc, 'sweep_m_rad', deg2rad(geom.htail_sweep_deg), 'fineness', NaN);
components(end+1) = struct('name',"fuselage",'type',"fuselage", ...
    'Swet', geom.lf_m * (geom.wf_max_m + geom.df_max_m) * 2, 'Lref', geom.lf_m, ...
    'tc', NaN, 'xmc', NaN, 'sweep_m_rad', NaN, ...
    'fineness', geom.lf_m / max(geom.wf_max_m, geom.df_max_m));

%% WEIGHT ITERATION
env = struct('rho', geom.rho, 'mu', geom.mu, 'a', geom.a);
dragParams = struct('Q', aero.Q, 'CD_misc', aero.CD_misc, 'CD_LP', aero.CD_LP);
fixedMasses = struct('payload_mass_kg', req.payload_mass_kg, ...
    'engine_mass_kg', mass.engine_mass_kg, 'avionics_mass_kg', mass.avionics_mass_kg, ...
    'landing_gear_mass_kg', mass.landing_gear_mass_kg, 'fuel_system_mass_kg', mass.fuel_system_mass_kg);

[W_final_N, S_final_m2, fuel_final_kg, tank_final_L, WeightHistory, converged, iter_converged, rel_error_final] = ...
    run_weight_iteration(req, assump, geom, fixedMasses, mass.W0_N, components, env, dragParams);

%% POST-CONVERGENCE PERFORMANCE ANALYSIS
components(1).Swet = S_final_m2 * 2;
components(1).Lref = sqrt(S_final_m2 / geom.ARw);

config = struct('W', W_final_N, 'Sref', S_final_m2, 'St', geom.Sh_m2, ...
    'ARw', geom.ARw, 'ARt', geom.ARh, 'components', components, ...
    'V_stall', req.V_stall_mps, 'V_cruise', req.V_design_mps, ...
    'fuel_mass_kg', fuel_final_kg, 'P_shaft_max', assump.P_shaft_max_W);

sweep_results = velocity_sweep(config, assump, env);

requirements_metadata = struct('Range_req_nmi_round_trip', req.Range_req_km * CONV.nmi_per_km, ...
    'Top_speed_req_mph', req.V_design_mps * CONV.mph_per_mps, ...
    'Stall_speed_req_mph', req.V_stall_mps * CONV.mph_per_mps, ...
    'Payload_weight_lb', req.payload_mass_kg * 2.20462262185);
sweep_results.Properties.UserData.requirements = requirements_metadata;


%% STABILITY ANALYSIS
cbar = sqrt(S_final_m2 / geom.ARw);
stabilityInputs = struct('cbar_m', cbar, 'x_ac_w_m', 0.25*cbar, ...
    'Sh_m2', geom.Sh_m2, 'eta_tail', stab.eta_tail, 'a_w', stab.a_w, 'a_t', stab.a_t, ...
    'deps_dalpha', stab.deps_dalpha, 'lt_m', stab.lt_frac*geom.lf_m);

compMassLoc = struct('name', {}, 'mass_kg', {}, 'x_m', {});
compMassLoc(end+1) = struct('name', "payload", 'mass_kg', fixedMasses.payload_mass_kg, 'x_m', cg.payload_x_m);
compMassLoc(end+1) = struct('name', "engine", 'mass_kg', fixedMasses.engine_mass_kg, 'x_m', cg.engine_x_m);
compMassLoc(end+1) = struct('name', "avionics", 'mass_kg', fixedMasses.avionics_mass_kg, 'x_m', cg.avionics_x_m);
compMassLoc(end+1) = struct('name', "fuel_system", 'mass_kg', fixedMasses.fuel_system_mass_kg, 'x_m', cg.fuel_system_x_m);
compMassLoc(end+1) = struct('name', "landing_gear", 'mass_kg', fixedMasses.landing_gear_mass_kg, 'x_m', cg.landing_gear_x_m);
compMassLoc(end+1) = struct('name', "fuel", 'mass_kg', fuel_final_kg, 'x_m', cg.fuel_x_m);

[x_cg_m, ~, CGTable] = compute_cg(compMassLoc);
x_np_m = compute_neutral_point(stabilityInputs, S_final_m2, geom.ARw);
[SM, ~] = compute_static_margin(x_cg_m, x_np_m, cbar, stab.SM_min, stab.SM_max);

%% POST-PROCESS RESULTS
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

V_80mph_mps = 80 / CONV.mph_per_mps;
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
results.xCG_m = x_cg_m;
results.xNP_m = x_np_m;
results.SM = SM;
results.CG_pctMAC = 100 * x_cg_m / cbar;
results.NP_pctMAC = 100 * x_np_m / cbar;
results.SM_pctMAC = 100 * SM;
results.CGTable = CGTable;
results.sweep_results = sweep_results;

%% PRINT RESULTS
printResults(results);

%% PLOTS
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
figure; hold on; grid on;
plot(V_vec, P_req, 'b-', 'LineWidth', 2, 'DisplayName', 'Power Required');
plot(V_vec, P_avail_vec, 'r-', 'LineWidth', 2, 'DisplayName', 'Power Available');
plot(V_vec, P_induced, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Induced Power');
plot(V_vec, P_parasite, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Parasite Power');

for k = 1:length(V_intersect_plot)
    Vmax = V_intersect_plot(k);
    Pmax = interp1(V_vec, P_req, Vmax);
    plot(Vmax, Pmax, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'y', 'HandleVisibility', 'off');
    text(Vmax, Pmax, sprintf('  V_{max} = %.1f m/s', Vmax), ...
        'VerticalAlignment', 'bottom', 'FontSize', 9);
end

plot(V_minP, Pmin, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'c', 'HandleVisibility', 'off');
text(V_minP, Pmin, sprintf('  Min Power = %.0f W', Pmin), ...
    'VerticalAlignment', 'top', 'FontSize', 9);

plot(V_Range_max, P_rng, 'kd', 'MarkerSize', 8, 'MarkerFaceColor', 'm', 'HandleVisibility', 'off');
text(V_Range_max, P_rng, sprintf('  V_{range} = %.1f m/s', V_Range_max), ...
    'VerticalAlignment', 'bottom', 'FontSize', 9);
xlabel('Velocity (m/s)');
ylabel('Power (W)');
title('Power Required and Power Available vs Velocity');
xlim([0.9*min(V_vec) 1.05*max(V_vec)]);
ylim([0 1.1*max(P_avail_vec)]);
legend('Location', 'best');
