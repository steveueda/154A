function printResults(results)

CONV.mph_per_mps = 2.236936;
CONV.nmi_per_km = 1 / 1.852;
percent_conversion = 100;
km_per_hr_per_mps = 3.6;

fprintf("\n");
fprintf("================================================================================\n");
fprintf("                    AIRCRAFT DESIGN ANALYSIS RESULTS\n");
fprintf("                    Run: %s\n", datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf("================================================================================\n");

% ============================================================================
% CONVERGENCE SUMMARY
% ============================================================================
fprintf("\n--- CONVERGENCE SUMMARY ---\n");
if results.converged
    fprintf("Converged at iteration %d\n", results.iter_converged);
    fprintf("Final relative error: %.3f%%\n", results.rel_error_final * 100);
else
    fprintf("Did not converge after %d iterations\n", results.max_iter);
    fprintf("Final relative error: %.3f%%\n", results.rel_error_final * 100);
end
fprintf("Final weight: %.2f N (%.2f kg)\n", results.W_final_N, results.W_final_kg);
fprintf("Final wing area: %.4f m^2\n", results.S_final_m2);

% ============================================================================
% WEIGHT BREAKDOWN
% ============================================================================
fprintf("\n--- WEIGHT BREAKDOWN ---\n");
fprintf("%-25s | %10s | %10s\n", "Component", "Mass (kg)", "Percent");
fprintf("%s\n", repmat('-', 1, 50));

if isfield(results, 'W_wing_kg') && ~isempty(results.W_wing_kg) && ~isnan(results.W_wing_kg)
    fprintf("STRUCTURAL COMPONENTS:\n");
    fprintf("%-25s | %10.3f | %9.1f%%\n", "  Wing", results.W_wing_kg, ...
        percent_conversion * results.W_wing_kg / results.W_final_kg);
    fprintf("%-25s | %10.3f | %9.1f%%\n", "  Fuselage", results.W_fuse_kg, ...
        percent_conversion * results.W_fuse_kg / results.W_final_kg);
    fprintf("%-25s | %10.3f | %9.1f%%\n", "  Horizontal Tail", results.W_ht_kg, ...
        percent_conversion * results.W_ht_kg / results.W_final_kg);
    fprintf("%-25s | %10.3f | %9.1f%%\n", "  Vertical Tail", results.W_vt_kg, ...
        percent_conversion * results.W_vt_kg / results.W_final_kg);
    fprintf("%-25s | %10.3f | %9.1f%%\n", "  Structural Subtotal", results.W_struct_kg, ...
        percent_conversion * results.W_struct_kg / results.W_final_kg);
    fprintf("  (Primary structure: wing + fuselage + tails)\n");
    fprintf("\n");
    fprintf("SECONDARY STRUCTURAL COMPONENTS:\n");
    fprintf("%-25s | %10.3f | %9.1f%%\n", "  Controls", results.W_controls_kg, ...
        percent_conversion * results.W_controls_kg / results.W_final_kg);
    fprintf("%-25s | %10.3f | %9.1f%%\n", "  Propulsion Install", results.W_prop_install_kg, ...
        percent_conversion * results.W_prop_install_kg / results.W_final_kg);
    if results.W_controls_kg > 0.5
        fprintf("  NOTE: Controls weight may be overestimated by Niccolai formula for small aircraft\n");
    end
    fprintf("\n");
end

CGTable = results.CGTable;
fprintf("FIXED COMPONENTS:\n");
for i = 1:height(CGTable)
    comp_name = CGTable.name(i);
    comp_mass = CGTable.mass_kg(i);
    if comp_mass > 0 && ~strcmpi(comp_name, "battery") && ~strcmpi(comp_name, "fuel")
        fprintf("%-25s | %10.3f | %9.1f%%\n", sprintf("  %s", comp_name), comp_mass, ...
            percent_conversion * comp_mass / results.W_final_kg);
    end
end

fprintf("\nENERGY STORAGE:\n");
for i = 1:height(CGTable)
    comp_name = CGTable.name(i);
    comp_mass = CGTable.mass_kg(i);
    if comp_mass > 0 && (strcmpi(comp_name, "battery") || strcmpi(comp_name, "fuel"))
        fprintf("%-25s | %10.3f | %9.1f%%\n", sprintf("  %s", comp_name), comp_mass, ...
            percent_conversion * comp_mass / results.W_final_kg);
    end
end

fprintf("%s\n", repmat('-', 1, 50));
fprintf("%-25s | %10.3f | %9.1f%%\n", "TOTAL GROSS WEIGHT", results.W_final_kg, 100.0);

% ============================================================================
% REQUIREMENTS VS PREDICTED PERFORMANCE
% ============================================================================
fprintf("\n--- REQUIREMENTS VS PREDICTED PERFORMANCE ---\n");
fprintf("%-25s | %-15s | %-15s | %-10s\n", "Requirement", "Predicted", "Margin", "Status");
fprintf("%s\n", repmat('-', 1, 75));

Range_req_nmi_round_trip = results.Range_req_km * CONV.nmi_per_km;

if isfield(results, 'Range_mission_nmi') && ~isnan(results.Range_mission_nmi)
    Range_mission_nmi = results.Range_mission_nmi;
    if Range_mission_nmi >= Range_req_nmi_round_trip
        Range_margin = Range_mission_nmi - Range_req_nmi_round_trip;
        Range_status = "PASS";
    else
        Range_margin = Range_mission_nmi - Range_req_nmi_round_trip;
        Range_status = "FAIL";
    end
    fprintf("%-25s | %-15s | %-15s | %-10s\n", ...
        sprintf("Range (round trip): %.1f nmi", Range_req_nmi_round_trip), ...
        sprintf("%.1f nmi", Range_mission_nmi), ...
        sprintf("%.1f nmi", Range_margin), Range_status);
else
    Range_one_leg_nmi = results.Range_best_km * CONV.nmi_per_km;
    Range_round_trip_nmi = Range_one_leg_nmi * 2;
    if Range_round_trip_nmi >= Range_req_nmi_round_trip
        Range_margin = Range_round_trip_nmi - Range_req_nmi_round_trip;
        Range_status = "PASS";
    else
        Range_margin = Range_round_trip_nmi - Range_req_nmi_round_trip;
        Range_status = "FAIL";
    end
    fprintf("%-25s | %-15s | %-15s | %-10s\n", ...
        sprintf("Range (round trip): %.1f nmi", Range_req_nmi_round_trip), ...
        sprintf("%.1f nmi", Range_round_trip_nmi), ...
        sprintf("%.1f nmi", Range_margin), Range_status);
    fprintf("   WARNING: Mission range not available, using theoretical best range speed\n");
end

V_top_req_mph = results.V_top_req_mph;
V_max_mph = results.V_max_mph;
if ~isnan(V_max_mph)
    if V_max_mph >= V_top_req_mph
        V_top_margin = V_max_mph - V_top_req_mph;
        V_top_status = "PASS";
    else
        V_top_margin = V_max_mph - V_top_req_mph;
        V_top_status = "FAIL";
    end
    fprintf("%-25s | %-15s | %-15s | %-10s\n", ...
        sprintf("Top speed: %.0f mph", V_top_req_mph), ...
        sprintf("%.1f mph", V_max_mph), ...
        sprintf("%.1f mph", V_top_margin), V_top_status);
else
    fprintf("%-25s | %-15s | %-15s | %-10s\n", ...
        sprintf("Top speed: %.0f mph", V_top_req_mph), ...
        "N/A", "N/A", "FAIL");
end

V_stall_req_mph = results.V_stall_req_mph;
V_stall_mph_computed = results.V_stall_mph;
if ~isnan(V_stall_mph_computed)
    if V_stall_mph_computed <= V_stall_req_mph
        V_stall_margin = V_stall_req_mph - V_stall_mph_computed;
        V_stall_status = "PASS";
    else
        V_stall_margin = V_stall_mph_computed - V_stall_req_mph;
        V_stall_status = "FAIL";
    end
    fprintf("%-25s | %-15s | %-15s | %-10s\n", ...
        sprintf("Stall speed: %.0f mph", V_stall_req_mph), ...
        sprintf("%.1f mph", V_stall_mph_computed), ...
        sprintf("%.1f mph", V_stall_margin), V_stall_status);
else
    fprintf("%-25s | %-15s | %-15s | %-10s\n", ...
        sprintf("Stall speed: %.0f mph", V_stall_req_mph), ...
        "N/A", "N/A", "FAIL");
end

% ============================================================================
% MISSION RANGE DETAILS
% ============================================================================
if isfield(results, 'Range_mission_nmi') && ~isnan(results.Range_mission_nmi)
    Range_one_leg_req_km = (results.Range_req_km / 2);
    Range_one_leg_req_nmi = Range_one_leg_req_km * CONV.nmi_per_km;
    Range_return_nmi = results.Range_mission_nmi - Range_one_leg_req_nmi;
    fprintf("\n--- MISSION RANGE DETAILS ---\n");
    fprintf("  Total mission range: %.1f nmi (round trip)\n", results.Range_mission_nmi);
    fprintf("  - Outbound leg: %.1f nmi at 80 mph\n", Range_one_leg_req_nmi);
    fprintf("  - Return leg: %.1f nmi at %.1f mph (best range speed)\n", Range_return_nmi, results.V_bestRange_return_mph);
end

% ============================================================================
% STABILITY SUMMARY
% ============================================================================
fprintf("\n--- STABILITY SUMMARY ---\n");
fprintf("Neutral point: %.3f m\n", results.xNP_m);

fprintf("\nOUTBOUND LEG (with payload, %.2f kg):\n", results.W_outbound_kg);
fprintf("  CG location: %.3f m\n", results.xCG_outbound_m);
fprintf("  Static margin: %.1f%% MAC\n", results.SM_outbound_pctMAC);

SM_min_target = 0.05;
SM_max_target = 0.15;
if results.SM_outbound > SM_min_target && results.SM_outbound < SM_max_target
    stability_verdict_outbound = "PASS";
elseif results.SM_outbound <= 0
    stability_verdict_outbound = "FAIL (insufficient margin)";
elseif results.SM_outbound >= SM_max_target
    stability_verdict_outbound = "PASS (but conservative)";
else
    stability_verdict_outbound = "FAIL (insufficient margin)";
end
fprintf("  Stability verdict: %s\n", stability_verdict_outbound);

fprintf("\nRETURN LEG (without payload, %.2f kg):\n", results.W_return_kg);
fprintf("  CG location: %.3f m\n", results.xCG_return_m);
fprintf("  Static margin: %.1f%% MAC\n", results.SM_return_pctMAC);

if results.SM_return > SM_min_target && results.SM_return < SM_max_target
    stability_verdict_return = "PASS";
elseif results.SM_return <= 0
    stability_verdict_return = "FAIL (insufficient margin)";
elseif results.SM_return >= SM_max_target
    stability_verdict_return = "PASS (but conservative)";
else
    stability_verdict_return = "FAIL (insufficient margin)";
end
fprintf("  Stability verdict: %s\n", stability_verdict_return);

% ============================================================================
% PERFORMANCE HIGHLIGHTS
% ============================================================================
fprintf("\n--- PERFORMANCE HIGHLIGHTS ---\n");

% Outbound leg
fprintf("OUTBOUND LEG (with payload):\n");
fprintf("  Top speed            : %6.1f m/s (%5.1f mph)\n", ...
    results.V_max_mps, results.V_max_mph);
fprintf("  Stall speed          : %6.1f m/s (%5.1f mph)\n", ...
    results.V_stall_mps, results.V_stall_mph);
fprintf("  Min power speed      : %6.1f m/s (%5.1f mph) | Preq = %6.1f W\n", ...
    results.V_minPower_mps, results.V_minPower_mph, results.Preq_min_W);

V_sweep = results.sweep_results.V_mps;
CL_sweep = results.sweep_results.CL;
CDtotal_sweep = results.sweep_results.CDtotal;
LD_sweep = CL_sweep ./ CDtotal_sweep;
[max_LD, idx_maxLD] = max(LD_sweep);
V_maxLD_mps = V_sweep(idx_maxLD);
V_maxLD_mph = V_maxLD_mps * CONV.mph_per_mps;

Range_one_leg_nmi_display = results.Range_best_km * CONV.nmi_per_km;
Range_round_trip_nmi_display = Range_one_leg_nmi_display * 2;
endurance_bestRange_hr = results.Range_best_km / (results.V_bestRange_mps * km_per_hr_per_mps);

fprintf("  Best range speed     : %6.1f m/s (%5.1f mph) | Range = %.1f km (%.1f nmi one leg, %.1f nmi round trip) | Endurance = %.2f hr\n", ...
    results.V_bestRange_mps, results.V_bestRange_mph, ...
    results.Range_best_km, Range_one_leg_nmi_display, Range_round_trip_nmi_display, endurance_bestRange_hr);
fprintf("  Max L/D              : %6.2f at %6.1f m/s (%5.1f mph)\n", ...
    max_LD, V_maxLD_mps, V_maxLD_mph);
fprintf("  Max ROC speed         : %6.1f m/s (%5.1f mph) | ROC = %.2f m/s, gamma = %.1f deg\n", ...
    results.V_maxROC_mps, results.V_maxROC_mph, ...
    results.ROC_max_mps, results.gamma_maxROC_deg);
if ~isnan(results.V_80mph_mps)
    fprintf("  At 80 mph requirement : %6.1f m/s (%5.1f mph) | Preq = %6.1f W, ExcessP = %6.1f W\n", ...
        results.V_80mph_mps, 80.0, ...
        results.Preq_80mph_W, results.ExcessP_80mph_W);
else
    fprintf("  At 80 mph requirement : N/A (outside sweep range)\n");
end

% Return leg
fprintf("\nRETURN LEG (without payload):\n");
fprintf("  Top speed            : %6.1f m/s (%5.1f mph)\n", ...
    results.V_max_return_mps, results.V_max_return_mph);
fprintf("  Stall speed          : %6.1f m/s (%5.1f mph)\n", ...
    results.V_stall_return_mps, results.V_stall_return_mph);
fprintf("  Min power speed      : %6.1f m/s (%5.1f mph) | Preq = %6.1f W\n", ...
    results.V_minPower_return_mps, results.V_minPower_return_mph, results.Preq_min_return_W);

Range_one_leg_return_nmi_display = results.Range_best_return_km * CONV.nmi_per_km;
Range_round_trip_return_nmi_display = Range_one_leg_return_nmi_display * 2;
endurance_bestRange_return_hr = results.Range_best_return_km / (results.V_bestRange_return_mps * km_per_hr_per_mps);

fprintf("  Best range speed     : %6.1f m/s (%5.1f mph) | Range = %.1f km (%.1f nmi one leg, %.1f nmi round trip) | Endurance = %.2f hr\n", ...
    results.V_bestRange_return_mps, results.V_bestRange_return_mph, ...
    results.Range_best_return_km, Range_one_leg_return_nmi_display, Range_round_trip_return_nmi_display, endurance_bestRange_return_hr);
fprintf("  Max L/D              : %6.2f at %6.1f m/s (%5.1f mph)\n", ...
    results.max_LD_return, results.V_maxLD_return_mps, results.V_maxLD_return_mph);
fprintf("  Max ROC speed         : %6.1f m/s (%5.1f mph) | ROC = %.2f m/s, gamma = %.1f deg\n", ...
    results.V_maxROC_return_mps, results.V_maxROC_return_mph, ...
    results.ROC_max_return_mps, results.gamma_maxROC_return_deg);
if ~isnan(results.Preq_80mph_return_W)
    fprintf("  At 80 mph requirement : %6.1f m/s (%5.1f mph) | Preq = %6.1f W, ExcessP = %6.1f W\n", ...
        results.V_80mph_mps, 80.0, ...
        results.Preq_80mph_return_W, results.ExcessP_80mph_return_W);
else
    fprintf("  At 80 mph requirement : N/A (outside sweep range)\n");
end

% ============================================================================
% SELECTED PERFORMANCE POINTS
% ============================================================================
fprintf("\n--- SELECTED PERFORMANCE POINTS ---\n");

% Outbound leg performance points
fprintf("OUTBOUND LEG (with payload):\n");
V_outbound = results.sweep_results.V_mps;
V_stall_outbound = results.V_stall_mps;
if ~isnan(results.V_max_mph)
    V_max_outbound = results.V_max_mph / CONV.mph_per_mps;
else
    V_max_outbound = NaN;
end
design_speed_mph = 80;
V_80mph_mps = design_speed_mph / CONV.mph_per_mps;
target_speeds_outbound = [V_stall_outbound, results.V_minPower_mps, ...
    results.V_bestRange_mps, results.V_maxROC_mps, V_80mph_mps, V_max_outbound];
target_labels_outbound = {"Stall speed", "Min power", "Best range", "Max ROC", "80 mph req", "Max speed"};

valid_idx_outbound = ~isnan(target_speeds_outbound);
target_speeds_outbound = target_speeds_outbound(valid_idx_outbound);
target_labels_outbound = target_labels_outbound(valid_idx_outbound);

Point_names_outbound = strings(length(target_speeds_outbound), 1);
V_mps_vec_outbound = zeros(length(target_speeds_outbound), 1);
V_mph_vec_outbound = zeros(length(target_speeds_outbound), 1);
CL_vec_outbound = zeros(length(target_speeds_outbound), 1);
CDp_vec_outbound = zeros(length(target_speeds_outbound), 1);
CDi_vec_outbound = zeros(length(target_speeds_outbound), 1);
CDtotal_vec_outbound = zeros(length(target_speeds_outbound), 1);
LD_vec_outbound = zeros(length(target_speeds_outbound), 1);
Drag_N_vec_outbound = zeros(length(target_speeds_outbound), 1);
Preq_W_vec_outbound = zeros(length(target_speeds_outbound), 1);
Pavail_W_vec_outbound = zeros(length(target_speeds_outbound), 1);
ExcessP_W_vec_outbound = zeros(length(target_speeds_outbound), 1);
ROC_mps_vec_outbound = zeros(length(target_speeds_outbound), 1);
gamma_deg_vec_outbound = zeros(length(target_speeds_outbound), 1);
Range_km_vec_outbound = zeros(length(target_speeds_outbound), 1);
Range_nmi_vec_outbound = zeros(length(target_speeds_outbound), 1);

for i = 1:length(target_speeds_outbound)
    V_targ = target_speeds_outbound(i);
    Point_names_outbound(i) = string(target_labels_outbound{i});
    
    if V_targ >= min(V_outbound) && V_targ <= max(V_outbound)
        V_mps_vec_outbound(i) = V_targ;
        V_mph_vec_outbound(i) = V_targ * CONV.mph_per_mps;
        CL_vec_outbound(i) = interp1(V_outbound, results.sweep_results.CL, V_targ);
        CDp_vec_outbound(i) = interp1(V_outbound, results.sweep_results.CDp, V_targ);
        CDi_vec_outbound(i) = interp1(V_outbound, results.sweep_results.CDi, V_targ);
        CDtotal_vec_outbound(i) = interp1(V_outbound, results.sweep_results.CDtotal, V_targ);
        LD_vec_outbound(i) = CL_vec_outbound(i) / CDtotal_vec_outbound(i);
        Drag_N_vec_outbound(i) = interp1(V_outbound, results.sweep_results.Drag_N, V_targ);
        Preq_W_vec_outbound(i) = interp1(V_outbound, results.sweep_results.Preq_W, V_targ);
        Pavail_W_vec_outbound(i) = results.sweep_results.Pavail_W(1);
        ExcessP_W_vec_outbound(i) = interp1(V_outbound, results.sweep_results.ExcessP_W, V_targ);
        ROC_mps_vec_outbound(i) = interp1(V_outbound, results.sweep_results.ROC_mps, V_targ);
        gamma_deg_vec_outbound(i) = interp1(V_outbound, results.sweep_results.gamma_deg, V_targ);
        Range_km_vec_outbound(i) = interp1(V_outbound, results.sweep_results.Range_km, V_targ);
        Range_nmi_vec_outbound(i) = Range_km_vec_outbound(i) * CONV.nmi_per_km;
    else
        [~, idx] = min(abs(V_outbound - V_targ));
        V_mps_vec_outbound(i) = V_outbound(idx);
        V_mph_vec_outbound(i) = V_mps_vec_outbound(i) * CONV.mph_per_mps;
        CL_vec_outbound(i) = results.sweep_results.CL(idx);
        CDp_vec_outbound(i) = results.sweep_results.CDp(idx);
        CDi_vec_outbound(i) = results.sweep_results.CDi(idx);
        CDtotal_vec_outbound(i) = results.sweep_results.CDtotal(idx);
        LD_vec_outbound(i) = CL_vec_outbound(i) / CDtotal_vec_outbound(i);
        Drag_N_vec_outbound(i) = results.sweep_results.Drag_N(idx);
        Preq_W_vec_outbound(i) = results.sweep_results.Preq_W(idx);
        Pavail_W_vec_outbound(i) = results.sweep_results.Pavail_W(idx);
        ExcessP_W_vec_outbound(i) = results.sweep_results.ExcessP_W(idx);
        ROC_mps_vec_outbound(i) = results.sweep_results.ROC_mps(idx);
        gamma_deg_vec_outbound(i) = results.sweep_results.gamma_deg(idx);
        Range_km_vec_outbound(i) = results.sweep_results.Range_km(idx);
        Range_nmi_vec_outbound(i) = Range_km_vec_outbound(i) * CONV.nmi_per_km;
    end
end

SelectedTable_outbound = table(Point_names_outbound, V_mps_vec_outbound, V_mph_vec_outbound, CL_vec_outbound, CDp_vec_outbound, CDi_vec_outbound, CDtotal_vec_outbound, LD_vec_outbound, ...
    Drag_N_vec_outbound, Preq_W_vec_outbound, Pavail_W_vec_outbound, ExcessP_W_vec_outbound, ROC_mps_vec_outbound, gamma_deg_vec_outbound, ...
    Range_km_vec_outbound, Range_nmi_vec_outbound, ...
    'VariableNames', {'Point', 'V_mps', 'V_mph', 'CL', 'CDp', 'CDi', 'CDtotal', 'L_D', ...
    'Drag_N', 'Preq_W', 'Pavail_W', 'ExcessP_W', 'ROC_mps', 'gamma_deg', 'Range_km', 'Range_nmi'});
disp(SelectedTable_outbound);

% Return leg performance points
fprintf("\nRETURN LEG (without payload):\n");
V_return = results.sweep_results_return.V_mps;
V_stall_return = results.V_stall_return_mps;
if ~isnan(results.V_max_return_mph)
    V_max_return = results.V_max_return_mph / CONV.mph_per_mps;
else
    V_max_return = NaN;
end
target_speeds_return = [V_stall_return, results.V_minPower_return_mps, ...
    results.V_bestRange_return_mps, results.V_maxROC_return_mps, V_80mph_mps, V_max_return];
target_labels_return = {"Stall speed", "Min power", "Best range", "Max ROC", "80 mph req", "Max speed"};

valid_idx_return = ~isnan(target_speeds_return);
target_speeds_return = target_speeds_return(valid_idx_return);
target_labels_return = target_labels_return(valid_idx_return);

Point_names_return = strings(length(target_speeds_return), 1);
V_mps_vec_return = zeros(length(target_speeds_return), 1);
V_mph_vec_return = zeros(length(target_speeds_return), 1);
CL_vec_return = zeros(length(target_speeds_return), 1);
CDp_vec_return = zeros(length(target_speeds_return), 1);
CDi_vec_return = zeros(length(target_speeds_return), 1);
CDtotal_vec_return = zeros(length(target_speeds_return), 1);
LD_vec_return = zeros(length(target_speeds_return), 1);
Drag_N_vec_return = zeros(length(target_speeds_return), 1);
Preq_W_vec_return = zeros(length(target_speeds_return), 1);
Pavail_W_vec_return = zeros(length(target_speeds_return), 1);
ExcessP_W_vec_return = zeros(length(target_speeds_return), 1);
ROC_mps_vec_return = zeros(length(target_speeds_return), 1);
gamma_deg_vec_return = zeros(length(target_speeds_return), 1);
Range_km_vec_return = zeros(length(target_speeds_return), 1);
Range_nmi_vec_return = zeros(length(target_speeds_return), 1);

for i = 1:length(target_speeds_return)
    V_targ = target_speeds_return(i);
    Point_names_return(i) = string(target_labels_return{i});
    
    if V_targ >= min(V_return) && V_targ <= max(V_return)
        V_mps_vec_return(i) = V_targ;
        V_mph_vec_return(i) = V_targ * CONV.mph_per_mps;
        CL_vec_return(i) = interp1(V_return, results.sweep_results_return.CL, V_targ);
        CDp_vec_return(i) = interp1(V_return, results.sweep_results_return.CDp, V_targ);
        CDi_vec_return(i) = interp1(V_return, results.sweep_results_return.CDi, V_targ);
        CDtotal_vec_return(i) = interp1(V_return, results.sweep_results_return.CDtotal, V_targ);
        LD_vec_return(i) = CL_vec_return(i) / CDtotal_vec_return(i);
        Drag_N_vec_return(i) = interp1(V_return, results.sweep_results_return.Drag_N, V_targ);
        Preq_W_vec_return(i) = interp1(V_return, results.sweep_results_return.Preq_W, V_targ);
        Pavail_W_vec_return(i) = results.sweep_results_return.Pavail_W(1);
        ExcessP_W_vec_return(i) = interp1(V_return, results.sweep_results_return.ExcessP_W, V_targ);
        ROC_mps_vec_return(i) = interp1(V_return, results.sweep_results_return.ROC_mps, V_targ);
        gamma_deg_vec_return(i) = interp1(V_return, results.sweep_results_return.gamma_deg, V_targ);
        Range_km_vec_return(i) = interp1(V_return, results.sweep_results_return.Range_km, V_targ);
        Range_nmi_vec_return(i) = Range_km_vec_return(i) * CONV.nmi_per_km;
    else
        [~, idx] = min(abs(V_return - V_targ));
        V_mps_vec_return(i) = V_return(idx);
        V_mph_vec_return(i) = V_mps_vec_return(i) * CONV.mph_per_mps;
        CL_vec_return(i) = results.sweep_results_return.CL(idx);
        CDp_vec_return(i) = results.sweep_results_return.CDp(idx);
        CDi_vec_return(i) = results.sweep_results_return.CDi(idx);
        CDtotal_vec_return(i) = results.sweep_results_return.CDtotal(idx);
        LD_vec_return(i) = CL_vec_return(i) / CDtotal_vec_return(i);
        Drag_N_vec_return(i) = results.sweep_results_return.Drag_N(idx);
        Preq_W_vec_return(i) = results.sweep_results_return.Preq_W(idx);
        Pavail_W_vec_return(i) = results.sweep_results_return.Pavail_W(idx);
        ExcessP_W_vec_return(i) = results.sweep_results_return.ExcessP_W(idx);
        ROC_mps_vec_return(i) = results.sweep_results_return.ROC_mps(idx);
        gamma_deg_vec_return(i) = results.sweep_results_return.gamma_deg(idx);
        Range_km_vec_return(i) = results.sweep_results_return.Range_km(idx);
        Range_nmi_vec_return(i) = Range_km_vec_return(i) * CONV.nmi_per_km;
    end
end

SelectedTable_return = table(Point_names_return, V_mps_vec_return, V_mph_vec_return, CL_vec_return, CDp_vec_return, CDi_vec_return, CDtotal_vec_return, LD_vec_return, ...
    Drag_N_vec_return, Preq_W_vec_return, Pavail_W_vec_return, ExcessP_W_vec_return, ROC_mps_vec_return, gamma_deg_vec_return, ...
    Range_km_vec_return, Range_nmi_vec_return, ...
    'VariableNames', {'Point', 'V_mps', 'V_mph', 'CL', 'CDp', 'CDi', 'CDtotal', 'L_D', ...
    'Drag_N', 'Preq_W', 'Pavail_W', 'ExcessP_W', 'ROC_mps', 'gamma_deg', 'Range_km', 'Range_nmi'});
disp(SelectedTable_return);

fprintf("\n================================================================================\n\n");

end
