function printResults(results)
CONV.mph_per_mps = 2.236936;
CONV.nmi_per_km = 1 / 1.852;

%% TITLE HEADER
fprintf("\n");
fprintf("================================================================================\n");
fprintf("                    AIRCRAFT DESIGN ANALYSIS RESULTS\n");
fprintf("                    Run: %s\n", datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf("================================================================================\n");

%% CONVERGENCE SUMMARY
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


%% REQUIREMENTS VS PREDICTED PERFORMANCE
fprintf("\n--- REQUIREMENTS VS PREDICTED PERFORMANCE ---\n");
fprintf("%-25s | %-15s | %-15s | %-10s\n", "Requirement", "Predicted", "Margin", "Status");
fprintf("%s\n", repmat('-', 1, 75));

% Range requirement (round trip)
Range_req_nmi_round_trip = results.Range_req_km * CONV.nmi_per_km;
Range_one_leg_nmi = results.Range_best_km * CONV.nmi_per_km;  % One leg (round trip = 2x)
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

% Top speed requirement
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

% Stall speed requirement (computed from CL = CLmax in sweep)
CL_max = 1.5;  % Maximum lift coefficient
V_stall_req_mph = results.V_stall_req_mph;
CL_check = results.sweep_results.CL;
V_check = results.sweep_results.V_mps;
idx_stall_check = find(CL_check >= CL_max, 1, 'last');
if isempty(idx_stall_check)
    [~, idx_stall_check] = min(abs(CL_check - CL_max));
    V_stall_mph_computed = V_check(idx_stall_check) * CONV.mph_per_mps;
elseif idx_stall_check == length(CL_check)
    V_stall_mph_computed = V_check(idx_stall_check) * CONV.mph_per_mps;
elseif idx_stall_check == 1
    V_stall_mph_computed = V_check(1) * CONV.mph_per_mps;
else
    V_stall_mps_computed = V_check(idx_stall_check) + (V_check(idx_stall_check+1) - V_check(idx_stall_check)) * ...
        (CL_max - CL_check(idx_stall_check)) / (CL_check(idx_stall_check+1) - CL_check(idx_stall_check));
    V_stall_mph_computed = V_stall_mps_computed * CONV.mph_per_mps;
end
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

%% KEY SPEEDS
fprintf("\n--- KEY SPEEDS ---\n");
fprintf("%-20s: %6.1f m/s (%5.1f mph) | Preq = %6.1f W\n", ...
    "Min power speed", results.V_minPower_mps, results.V_minPower_mph, results.Preq_min_W);
Range_one_leg_nmi_display = results.Range_best_km * CONV.nmi_per_km;
Range_round_trip_nmi_display = Range_one_leg_nmi_display * 2;
fprintf("%-20s: %6.1f m/s (%5.1f mph) | Range = %.1f km (%.1f nmi one leg, %.1f nmi round trip)\n", ...
    "Best range speed", results.V_bestRange_mps, results.V_bestRange_mph, ...
    results.Range_best_km, Range_one_leg_nmi_display, Range_round_trip_nmi_display);
fprintf("%-20s: %6.1f m/s (%5.1f mph) | ROC = %.2f m/s, gamma = %.1f deg\n", ...
    "Max ROC speed", results.V_maxROC_mps, results.V_maxROC_mph, ...
    results.ROC_max_mps, results.gamma_maxROC_deg);
if ~isnan(results.V_80mph_mps)
    fprintf("%-20s: %6.1f m/s (%5.1f mph) | Preq = %6.1f W, ExcessP = %6.1f W\n", ...
        "At 80 mph requirement", results.V_80mph_mps, 80.0, ...
        results.Preq_80mph_W, results.ExcessP_80mph_W);
else
    fprintf("%-20s: %6s (%5s) | 80 mph outside sweep range\n", ...
        "At 80 mph requirement", "N/A", "N/A");
end

%% STABILITY SUMMARY
fprintf("\n--- STABILITY SUMMARY ---\n");
fprintf("CG location: %.3f m (%.1f%% MAC)\n", results.xCG_m, results.CG_pctMAC);
fprintf("Neutral point: %.3f m (%.1f%% MAC)\n", results.xNP_m, results.NP_pctMAC);
fprintf("Static margin: %.1f%% MAC\n", results.SM_pctMAC);

% Stability verdict
if results.SM > 0.05 && results.SM < 0.15
    stability_verdict = "PASS";
elseif results.SM <= 0
    stability_verdict = "FAIL (insufficient margin)";
elseif results.SM >= 0.15
    stability_verdict = "PASS (but conservative)";
else
    stability_verdict = "FAIL (insufficient margin)";
end
fprintf("Stability verdict: %s\n", stability_verdict);


%% COMPACT VELOCITY SWEEP TABLE
fprintf("\n--- SELECTED PERFORMANCE POINTS ---\n");

% Find indices for key speeds
V = results.sweep_results.V_mps;
CL = results.sweep_results.CL;
CL_max = 1.5;  % Maximum lift coefficient

% Compute stall speed from CL = CLmax
% Find where CL crosses CLmax (CL decreases with increasing V)
idx_stall = find(CL >= CL_max, 1, 'last');  % Last point where CL >= CLmax
if isempty(idx_stall)
    % CL never reaches CLmax - find closest point
    [~, idx_stall] = min(abs(CL - CL_max));
    V_stall_mps = V(idx_stall);
elseif idx_stall == length(CL)
    % CL is always >= CLmax at end - use last point
    V_stall_mps = V(idx_stall);
elseif idx_stall == 1
    % CL is always >= CLmax - use first point
    V_stall_mps = V(1);
else
    % Interpolate between idx_stall and idx_stall+1
    V_stall_mps = V(idx_stall) + (V(idx_stall+1) - V(idx_stall)) * ...
        (CL_max - CL(idx_stall)) / (CL(idx_stall+1) - CL(idx_stall));
end

if ~isnan(results.V_max_mph)
    V_max_mps = results.V_max_mph / CONV.mph_per_mps;
else
    V_max_mps = NaN;
end
V_80mph_mps = 80 / CONV.mph_per_mps;  % Required top speed (80 mph in m/s)
target_speeds = [V_stall_mps, results.V_minPower_mps, ...
    results.V_bestRange_mps, results.V_maxROC_mps, V_80mph_mps, V_max_mps];
target_labels = {"Stall speed", "Min power", "Best range", "Max ROC", "80 mph req", "Max speed"};

% Remove NaN targets
valid_idx = ~isnan(target_speeds);
target_speeds = target_speeds(valid_idx);
target_labels = target_labels(valid_idx);

% Interpolate or find nearest for each target
Point_names = strings(length(target_speeds), 1);
V_mps_vec = zeros(length(target_speeds), 1);
V_mph_vec = zeros(length(target_speeds), 1);
CL_vec = zeros(length(target_speeds), 1);
CDp_vec = zeros(length(target_speeds), 1);
CDi_vec = zeros(length(target_speeds), 1);
CDtotal_vec = zeros(length(target_speeds), 1);
Drag_N_vec = zeros(length(target_speeds), 1);
Preq_W_vec = zeros(length(target_speeds), 1);
Pavail_W_vec = zeros(length(target_speeds), 1);
ExcessP_W_vec = zeros(length(target_speeds), 1);
ROC_mps_vec = zeros(length(target_speeds), 1);
gamma_deg_vec = zeros(length(target_speeds), 1);
Range_km_vec = zeros(length(target_speeds), 1);
Range_nmi_vec = zeros(length(target_speeds), 1);

for i = 1:length(target_speeds)
    V_targ = target_speeds(i);
    Point_names(i) = string(target_labels{i});
    
    % Special handling for stall speed - use CL = CLmax
    if strcmp(target_labels{i}, "Stall speed")
        V_mps_vec(i) = V_stall_mps;
        V_mph_vec(i) = V_stall_mps * CONV.mph_per_mps;
        CL_vec(i) = CL_max;  % Set explicitly to CLmax
        % Interpolate other quantities at stall speed
        if V_stall_mps >= min(V) && V_stall_mps <= max(V)
            CDp_vec(i) = interp1(V, results.sweep_results.CDp, V_stall_mps);
            CDi_vec(i) = interp1(V, results.sweep_results.CDi, V_stall_mps);
            CDtotal_vec(i) = interp1(V, results.sweep_results.CDtotal, V_stall_mps);
            Drag_N_vec(i) = interp1(V, results.sweep_results.Drag_N, V_stall_mps);
            Preq_W_vec(i) = interp1(V, results.sweep_results.Preq_W, V_stall_mps);
            Pavail_W_vec(i) = results.sweep_results.Pavail_W(1);
            ExcessP_W_vec(i) = interp1(V, results.sweep_results.ExcessP_W, V_stall_mps);
            ROC_mps_vec(i) = interp1(V, results.sweep_results.ROC_mps, V_stall_mps);
            gamma_deg_vec(i) = interp1(V, results.sweep_results.gamma_deg, V_stall_mps);
            Range_km_vec(i) = interp1(V, results.sweep_results.Range_km, V_stall_mps);
        else
            % Find nearest if outside range
            [~, idx] = min(abs(V - V_stall_mps));
            CDp_vec(i) = results.sweep_results.CDp(idx);
            CDi_vec(i) = results.sweep_results.CDi(idx);
            CDtotal_vec(i) = results.sweep_results.CDtotal(idx);
            Drag_N_vec(i) = results.sweep_results.Drag_N(idx);
            Preq_W_vec(i) = results.sweep_results.Preq_W(idx);
            Pavail_W_vec(i) = results.sweep_results.Pavail_W(idx);
            ExcessP_W_vec(i) = results.sweep_results.ExcessP_W(idx);
            ROC_mps_vec(i) = results.sweep_results.ROC_mps(idx);
            gamma_deg_vec(i) = results.sweep_results.gamma_deg(idx);
            Range_km_vec(i) = results.sweep_results.Range_km(idx);
        end
        Range_nmi_vec(i) = Range_km_vec(i) * CONV.nmi_per_km;
    elseif V_targ >= min(V) && V_targ <= max(V)
        % Interpolate for other points
        V_mps_vec(i) = V_targ;
        V_mph_vec(i) = V_targ * CONV.mph_per_mps;
        CL_vec(i) = interp1(V, results.sweep_results.CL, V_targ);
        CDp_vec(i) = interp1(V, results.sweep_results.CDp, V_targ);
        CDi_vec(i) = interp1(V, results.sweep_results.CDi, V_targ);
        CDtotal_vec(i) = interp1(V, results.sweep_results.CDtotal, V_targ);
        Drag_N_vec(i) = interp1(V, results.sweep_results.Drag_N, V_targ);
        Preq_W_vec(i) = interp1(V, results.sweep_results.Preq_W, V_targ);
        Pavail_W_vec(i) = results.sweep_results.Pavail_W(1);
        ExcessP_W_vec(i) = interp1(V, results.sweep_results.ExcessP_W, V_targ);
        ROC_mps_vec(i) = interp1(V, results.sweep_results.ROC_mps, V_targ);
        gamma_deg_vec(i) = interp1(V, results.sweep_results.gamma_deg, V_targ);
        Range_km_vec(i) = interp1(V, results.sweep_results.Range_km, V_targ);
        Range_nmi_vec(i) = Range_km_vec(i) * CONV.nmi_per_km;
    else
        % Find nearest
        [~, idx] = min(abs(V - V_targ));
        V_mps_vec(i) = V(idx);
        V_mph_vec(i) = V_mps_vec(i) * CONV.mph_per_mps;
        CL_vec(i) = results.sweep_results.CL(idx);
        CDp_vec(i) = results.sweep_results.CDp(idx);
        CDi_vec(i) = results.sweep_results.CDi(idx);
        CDtotal_vec(i) = results.sweep_results.CDtotal(idx);
        Drag_N_vec(i) = results.sweep_results.Drag_N(idx);
        Preq_W_vec(i) = results.sweep_results.Preq_W(idx);
        Pavail_W_vec(i) = results.sweep_results.Pavail_W(idx);
        ExcessP_W_vec(i) = results.sweep_results.ExcessP_W(idx);
        ROC_mps_vec(i) = results.sweep_results.ROC_mps(idx);
        gamma_deg_vec(i) = results.sweep_results.gamma_deg(idx);
        Range_km_vec(i) = results.sweep_results.Range_km(idx);
        Range_nmi_vec(i) = Range_km_vec(i) * CONV.nmi_per_km;
    end
end

% Create and display table
SelectedTable = table(Point_names, V_mps_vec, V_mph_vec, CL_vec, CDp_vec, CDi_vec, CDtotal_vec, ...
    Drag_N_vec, Preq_W_vec, Pavail_W_vec, ExcessP_W_vec, ROC_mps_vec, gamma_deg_vec, ...
    Range_km_vec, Range_nmi_vec, ...
    'VariableNames', {'Point', 'V_mps', 'V_mph', 'CL', 'CDp', 'CDi', 'CDtotal', ...
    'Drag_N', 'Preq_W', 'Pavail_W', 'ExcessP_W', 'ROC_mps', 'gamma_deg', 'Range_km', 'Range_nmi'});
disp(SelectedTable);



fprintf("\n================================================================================\n\n");

end

