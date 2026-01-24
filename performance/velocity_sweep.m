function PerfTable = velocity_sweep(config, assump, env)
% Post-convergence performance analysis. Does not modify weight, geometry, or CG.

% Velocity sweep bounds (adaptive to capture full performance envelope)
V_stall = config.V_stall;
V_guess_cruise = config.V_cruise;
V_min = 0.6 * V_stall;      % Start below stall for complete curve
V_max = 1.5 * V_guess_cruise;  % Extend beyond cruise to find max speed
N_points = 150;             % High resolution for smooth curves

% Initial sweep
V = linspace(V_min, V_max, N_points)';
[PerfTable, ExcessP] = compute_sweep(V, config, assump, env);

% Check for power intersections (need 2 for max speed determination)
idx_cross = find(diff(sign(ExcessP)) ~= 0);
if length(idx_cross) < 2
    % Expand V_max and recompute once to capture both intersections
    V_max = 2.0 * V_guess_cruise;
    V = linspace(V_min, V_max, N_points)';
    [PerfTable, ExcessP] = compute_sweep(V, config, assump, env);
end

% Print performance highlights (requirements passed via PerfTable metadata)
printPerformanceHighlights(PerfTable, config, assump);

end

function [PerfTable, ExcessP] = compute_sweep(V, config, assump, env)
% Helper function to compute performance sweep

q = 0.5 * env.rho .* V.^2;
CL = config.W ./ (q .* config.Sref);

% Parasite drag
dragParams = struct('Q', assump.Q, 'CD_misc', assump.CD_misc, 'CD_LP', assump.CD_LP);
aero = struct('ARw', config.ARw, 'ew', assump.ew);
CDp = zeros(size(V));
for j = 1:length(V)
    drag = compute_drag(V(j), config.W, config.Sref, config.components, env, aero, dragParams);
    CDp(j) = drag.CDp;
end

% Induced drag (wing only)
CDi = CL.^2 ./ (pi * config.ARw * assump.ew);
CDtotal = CDp + CDi;
Drag = q .* config.Sref .* CDtotal;

% Power
Preq = Drag .* V;
Pavail = assump.eta_prop * config.P_shaft_max;
ExcessP = Pavail - Preq;

% Climb
ROC = ExcessP ./ config.W;
ROC_over_V = ROC ./ V;
ROC_over_V(ROC_over_V > 1) = 1;
ROC_over_V(ROC_over_V < -1) = -1;
gamma_deg = rad2deg(asin(ROC_over_V));

% Range (BSFC model)
Pshaft_kW = (Preq ./ assump.eta_prop) / 1000;
mdot_fuel_kg_per_hr = (assump.BSFC_g_per_kWh / 1000) .* Pshaft_kW;
endurance_hr = config.fuel_mass_kg ./ mdot_fuel_kg_per_hr;
range_km = V .* endurance_hr * 3.6;

PerfTable = table(V, q, CL, CDp, CDi, CDtotal, Drag, Preq, ...
    Pavail*ones(size(V)), ExcessP, ROC, gamma_deg, range_km, ...
    'VariableNames', {'V_mps','q_Pa','CL','CDp','CDi','CDtotal','Drag_N','Preq_W', ...
    'Pavail_W','ExcessP_W','ROC_mps','gamma_deg','Range_km'});

end

function printPerformanceHighlights(PerfTable, config, assump)
% Prints formatted performance highlights and requirement checks

% Unit conversion constants
CONV = struct();
CONV.mph_per_mps = 2.236936;
CONV.nmi_per_km = 1 / 1.852;
CONV.ft_per_m = 3.28084;
CONV.fpm_per_mps = 196.8504;

% Extract requirements from PerfTable metadata (set in main script)
if isfield(PerfTable.Properties, 'UserData') && isfield(PerfTable.Properties.UserData, 'requirements')
    requirements = PerfTable.Properties.UserData.requirements;
else
    % Fallback defaults (should not occur in normal operation)
    requirements = struct('Range_req_nmi_round_trip', 30, 'Top_speed_req_mph', 80, ...
        'Stall_speed_req_mph', 40, 'Payload_weight_lb', 4.5);
end

CL_max = assump.CL_max;

V = PerfTable.V_mps;
Preq = PerfTable.Preq_W;
Pavail = PerfTable.Pavail_W(1);
ExcessP = PerfTable.ExcessP_W;
ROC = PerfTable.ROC_mps;
gamma_deg = PerfTable.gamma_deg;
Range_km = PerfTable.Range_km;
CL = PerfTable.CL;

% 1) Stall speed estimate
if all(CL < CL_max)
    Vstall_mps = NaN;
    Vstall_mph = NaN;
    stall_status = " (stall not captured, increase low end of sweep)";
elseif all(CL >= CL_max)
    Vstall_mps = NaN;
    Vstall_mph = NaN;
    stall_status = " (all points stalled, increase speed range)";
else
    idx_stall = find(CL >= CL_max, 1);
    if idx_stall == 1
        Vstall_mps = V(1);
    else
        % Linear interpolation
        Vstall_mps = V(idx_stall-1) + (V(idx_stall) - V(idx_stall-1)) * ...
            (CL_max - CL(idx_stall-1)) / (CL(idx_stall) - CL(idx_stall-1));
    end
    Vstall_mph = Vstall_mps * CONV.mph_per_mps;
    stall_status = "";
end

% 2) Max speed (power intersections)
idx_cross = find(diff(sign(ExcessP)) ~= 0);
V_intersect_list = [];
for i = 1:length(idx_cross)
    idx = idx_cross(i);
    V_int = V(idx) + (V(idx+1) - V(idx)) * ...
        (-ExcessP(idx)) / (ExcessP(idx+1) - ExcessP(idx));
    V_intersect_list(end+1) = V_int;
end
if length(V_intersect_list) >= 2
    Vmax_mps = max(V_intersect_list);
    Vmax_mph = Vmax_mps * CONV.mph_per_mps;
    max_speed_warning = "";
else
    if ~isempty(V_intersect_list)
        Vmax_mps = max(V_intersect_list);
        Vmax_mph = Vmax_mps * CONV.mph_per_mps;
        max_speed_warning = " (WARNING: fewer than 2 intersections)";
    else
        Vmax_mps = NaN;
        Vmax_mph = NaN;
        max_speed_warning = " (WARNING: no intersections found)";
    end
end

% 3) Minimum power speed
[Pmin, idx_minP] = min(Preq);
V_minPower_mps = V(idx_minP);
V_minPower_mph = V_minPower_mps * CONV.mph_per_mps;

% 4) Best range speed (maximize V/Preq metric)
range_metric = V ./ Preq;
[~, idx_bestRange] = max(range_metric);
V_bestRange_mps = V(idx_bestRange);
V_bestRange_mph = V_bestRange_mps * CONV.mph_per_mps;
Range_one_leg_nmi = Range_km(idx_bestRange) * CONV.nmi_per_km;
Range_round_trip_nmi = Range_one_leg_nmi * 2;  % Round trip = 2x one leg

% 5) Best climb (maximum rate of climb)
[ROC_max_mps, idx_bestROC] = max(ROC);
ROC_max_fpm = ROC_max_mps * CONV.fpm_per_mps;
V_at_ROCmax_mph = V(idx_bestROC) * CONV.mph_per_mps;
gamma_deg_at_ROCmax = gamma_deg(idx_bestROC);

% 6) Power margin at required top speed
V_req_top_mps = requirements.Top_speed_req_mph / CONV.mph_per_mps;
if V_req_top_mps >= min(V) && V_req_top_mps <= max(V)
    Preq_at_80mph = interp1(V, Preq, V_req_top_mps);
    ExcessP_at_80mph = interp1(V, ExcessP, V_req_top_mps);
else
    Preq_at_80mph = NaN;
    ExcessP_at_80mph = NaN;
end

% 7) Range requirement check (round trip requirement)
Range_req_km_round_trip = requirements.Range_req_nmi_round_trip / CONV.nmi_per_km;
Range_one_leg_req_km = Range_req_km_round_trip / 2;  % Split round trip in half
Range_one_leg_km = Range_one_leg_nmi / CONV.nmi_per_km;
if Range_one_leg_km >= Range_one_leg_req_km
    range_margin_pct = ((Range_one_leg_km / Range_one_leg_req_km) - 1) * 100;
    range_check = sprintf("PASS (+%.1f%%)", range_margin_pct);
else
    range_margin_pct = (1 - (Range_one_leg_km / Range_one_leg_req_km)) * 100;
    range_check = sprintf("FAIL (-%.1f%%)", range_margin_pct);
end

% 8) Top speed requirement check
if ~isnan(Vmax_mph)
    if Vmax_mph >= requirements.Top_speed_req_mph
        top_speed_margin = Vmax_mph - requirements.Top_speed_req_mph;
        top_speed_check = sprintf("PASS (+%.1f mph)", top_speed_margin);
    else
        top_speed_margin = requirements.Top_speed_req_mph - Vmax_mph;
        top_speed_check = sprintf("FAIL (-%.1f mph)", top_speed_margin);
    end
else
    top_speed_check = "FAIL (no max speed found)";
end

% 9) Stall speed requirement check
if ~isnan(Vstall_mph)
    if Vstall_mph <= requirements.Stall_speed_req_mph
        stall_margin = requirements.Stall_speed_req_mph - Vstall_mph;
        stall_check = sprintf("PASS (+%.1f mph margin)", stall_margin);
    else
        stall_margin = Vstall_mph - requirements.Stall_speed_req_mph;
        stall_check = sprintf("FAIL (+%.1f mph over)", stall_margin);
    end
else
    stall_check = sprintf("FAIL (%s)", stall_status);
end

% Print formatted output
fprintf("\n==================== PERFORMANCE HIGHLIGHTS ====================\n");
fprintf("A) Requirements:\n");
fprintf("   Range (round trip): %.0f nmi\n", requirements.Range_req_nmi_round_trip);
fprintf("   Top speed: %.0f mph\n", requirements.Top_speed_req_mph);
fprintf("   Stall speed: %.0f mph\n", requirements.Stall_speed_req_mph);
fprintf("   Payload: %.1f lb\n", requirements.Payload_weight_lb);

fprintf("\nB) Predicted Performance:\n");
if ~isnan(Vstall_mph)
    fprintf("   Stall speed: %.1f mph (%.1f m/s)%s\n", Vstall_mph, Vstall_mps, stall_status);
else
    fprintf("   Stall speed: N/A%s\n", stall_status);
end
if ~isnan(Vmax_mph)
    fprintf("   Max speed: %.1f mph (%.1f m/s)%s\n", Vmax_mph, Vmax_mps, max_speed_warning);
else
    fprintf("   Max speed: N/A%s\n", max_speed_warning);
end
fprintf("   Min power speed: %.1f mph (%.1f m/s), Preq = %.0f W\n", ...
    V_minPower_mph, V_minPower_mps, Pmin);
fprintf("   Best range speed: %.1f mph (%.1f m/s)\n", V_bestRange_mph, V_bestRange_mps);
fprintf("   Range at best speed: %.1f nmi (one leg), %.1f nmi (round trip)\n", ...
    Range_one_leg_nmi, Range_round_trip_nmi);
fprintf("   Best climb: %.1f m/s (%.0f fpm) at %.1f mph, %.1f deg\n", ...
    ROC_max_mps, ROC_max_fpm, V_at_ROCmax_mph, gamma_deg_at_ROCmax);
if ~isnan(Preq_at_80mph)
    fprintf("   Power at 80 mph: Preq = %.0f W, ExcessP = %.0f W\n", ...
        Preq_at_80mph, ExcessP_at_80mph);
else
    fprintf("   Power at 80 mph: 80 mph outside sweep range\n");
end

fprintf("\nC) Requirement Checks:\n");
fprintf("   Range: %s\n", range_check);
fprintf("   Top speed: %s\n", top_speed_check);
fprintf("   Stall speed: %s\n", stall_check);
fprintf("============================================================\n\n");

end
