function PerfTable = velocity_sweep(config, assump, env)

V_stall = config.V_stall;
V_guess_cruise = config.V_cruise;
V_min_fraction = assump.V_min_fraction;
V_max_fraction = assump.V_max_fraction;
V_max_expanded_fraction = assump.V_max_fraction;
N_sweep_points = 150;
min_intersections_required = 2;
tolerance_nmi = 0.01;

V_min = V_min_fraction * V_stall;
V_max = V_max_fraction * V_guess_cruise;
V = linspace(V_min, V_max, N_sweep_points)';
[PerfTable, ExcessP] = compute_sweep(V, config, assump, env);

idx_cross = find(diff(sign(ExcessP)) ~= 0);
if length(idx_cross) < min_intersections_required
    V_max = V_max_expanded_fraction * V_guess_cruise;
    V = linspace(V_min, V_max, N_sweep_points)';
    [PerfTable, ExcessP] = compute_sweep(V, config, assump, env);
end

Range_mission_nmi = printPerformanceHighlights(PerfTable, config, assump, env);

if isempty(PerfTable.Properties.UserData) || ~isstruct(PerfTable.Properties.UserData)
    PerfTable.Properties.UserData = struct('mission_range_nmi', Range_mission_nmi);
else
    PerfTable.Properties.UserData.mission_range_nmi = Range_mission_nmi;
end

if isnan(Range_mission_nmi)
    warning('Mission range calculation returned NaN - check velocity sweep and energy calculations');
elseif isfield(PerfTable.Properties.UserData, 'mission_range_nmi')
    stored_range = PerfTable.Properties.UserData.mission_range_nmi;
    if abs(stored_range - Range_mission_nmi) > tolerance_nmi
        warning('Mission range stored value (%.2f) differs from calculated (%.2f)', ...
            stored_range, Range_mission_nmi);
    end
end

end

function [PerfTable, ExcessP] = compute_sweep(V, config, assump, env)

q = 0.5 * env.rho .* V.^2;
CL = config.W ./ (q .* config.Sref);

dragParams = struct('Q', assump.Q, 'CD_misc', assump.CD_misc, 'CD_LP', assump.CD_LP);
aero = struct('ARw', config.ARw, 'ew', assump.ew);
CDp = zeros(size(V));
for j = 1:length(V)
    drag = compute_drag(V(j), config.W, config.Sref, config.components, env, aero, dragParams);
    CDp(j) = drag.CDp;
end

CDi = CL.^2 ./ (pi * config.ARw * assump.ew);
CDtotal = CDp + CDi;
Drag = q .* config.Sref .* CDtotal;

Preq = Drag .* V;
Pavail = assump.eta_prop * config.P_shaft_max;
ExcessP = Pavail - Preq;

ROC = ExcessP ./ config.W;
ROC_over_V = ROC ./ V;
ROC_over_V = max(-1, min(1, ROC_over_V));
gamma_deg = rad2deg(asin(ROC_over_V));

km_per_hr_per_mps = 3.6;
grams_per_kg = 1000;
watts_per_kW = 1000;

if strcmpi(assump.prop_type, 'electric')
    P_electric_W = Preq ./ (assump.eta_prop * assump.eta_electric);
    E_available_Wh = config.fuel_mass_kg * assump.battery_energy_density_Wh_per_kg;
    endurance_hr = E_available_Wh ./ P_electric_W;
    range_km = V .* endurance_hr * km_per_hr_per_mps;
else
    Pshaft_kW = (Preq ./ assump.eta_prop) / watts_per_kW;
    BSFC_kg_per_kWh = assump.BSFC_g_per_kWh / grams_per_kg;
    mdot_fuel_kg_per_hr = BSFC_kg_per_kWh .* Pshaft_kW;
    endurance_hr = config.fuel_mass_kg ./ mdot_fuel_kg_per_hr;
    range_km = V .* endurance_hr * km_per_hr_per_mps;
end

PerfTable = table(V, q, CL, CDp, CDi, CDtotal, Drag, Preq, ...
    Pavail*ones(size(V)), ExcessP, ROC, gamma_deg, range_km, ...
    'VariableNames', {'V_mps','q_Pa','CL','CDp','CDi','CDtotal','Drag_N','Preq_W', ...
    'Pavail_W','ExcessP_W','ROC_mps','gamma_deg','Range_km'});

end

function Range_mission_nmi = printPerformanceHighlights(PerfTable, config, assump, env)

CONV = struct();
CONV.mph_per_mps = 2.236936;
CONV.nmi_per_km = 1 / 1.852;
CONV.ft_per_m = 3.28084;
CONV.fpm_per_mps = 196.8504;
CONV.meters_per_km = 1000;
CONV.seconds_per_hour = 3600;
CONV.km_per_hr_per_mps = 3.6;
CONV.grams_per_kg = 1000;
CONV.watts_per_kW = 1000;

if isfield(PerfTable.Properties, 'UserData') && isfield(PerfTable.Properties.UserData, 'requirements')
    requirements = PerfTable.Properties.UserData.requirements;
else
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
        Vstall_mps = V(idx_stall-1) + (V(idx_stall) - V(idx_stall-1)) * ...
            (CL_max - CL(idx_stall-1)) / (CL(idx_stall) - CL(idx_stall-1));
    end
    Vstall_mph = Vstall_mps * CONV.mph_per_mps;
    stall_status = "";
end

idx_cross = find(diff(sign(ExcessP)) ~= 0);
V_intersect_list = [];
for i = 1:length(idx_cross)
    idx = idx_cross(i);
    V_int = V(idx) + (V(idx+1) - V(idx)) * ...
        (-ExcessP(idx)) / (ExcessP(idx+1) - ExcessP(idx));
    V_intersect_list(end+1) = V_int;
end
min_intersections_required = 2;
if length(V_intersect_list) >= min_intersections_required
    Vmax_mps = max(V_intersect_list);
    Vmax_mph = Vmax_mps * CONV.mph_per_mps;
    max_speed_warning = "";
elseif ~isempty(V_intersect_list)
    Vmax_mps = max(V_intersect_list);
    Vmax_mph = Vmax_mps * CONV.mph_per_mps;
    max_speed_warning = " (WARNING: fewer than 2 intersections)";
else
    Vmax_mps = NaN;
    Vmax_mph = NaN;
    max_speed_warning = " (WARNING: no intersections found)";
end

[Pmin, idx_minP] = min(Preq);
V_minPower_mps = V(idx_minP);
V_minPower_mph = V_minPower_mps * CONV.mph_per_mps;

range_metric = V ./ Preq;
[~, idx_bestRange] = max(range_metric);
V_bestRange_mps = V(idx_bestRange);
V_bestRange_mph = V_bestRange_mps * CONV.mph_per_mps;
Range_one_leg_nmi = Range_km(idx_bestRange) * CONV.nmi_per_km;
Range_round_trip_nmi = Range_one_leg_nmi * 2;

Range_mission_nmi = NaN;
Range_return_nmi = NaN;
V_req_top_mps = requirements.Top_speed_req_mph / CONV.mph_per_mps;
V_outbound_mps = V_req_top_mps;

if V_outbound_mps >= min(V) && V_outbound_mps <= max(V) && ~isnan(Preq(idx_bestRange))
    Preq_outbound = interp1(V, Preq, V_outbound_mps);
    Range_one_leg_req_km = (requirements.Range_req_nmi_round_trip / 2) / CONV.nmi_per_km;
    t_outbound_hr = (Range_one_leg_req_km * CONV.meters_per_km) / (V_outbound_mps * CONV.seconds_per_hour);
    
    if isfield(config, 'W_return') && isfield(config, 'payload_mass_kg')
        W_return_N = config.W_return;
        payload_weight_N = config.payload_mass_kg * 9.80665;
        
        dragParams = struct('Q', assump.Q, 'CD_misc', assump.CD_misc, 'CD_LP', assump.CD_LP);
        aero = struct('ARw', config.ARw, 'ew', assump.ew);
        
        V_stall_return_est = sqrt(2 * W_return_N / (env.rho * config.Sref * assump.CL_max));
        V_min_search_return = max(V_stall_return_est * 1.1, 10);
        V_max_search_return = V_outbound_mps * 0.95;
        [V_bestRange_return_mps, Preq_return] = find_best_range_speed(W_return_N, config.Sref, config.components, env, aero, dragParams, assump, V_min_search_return, V_max_search_return);
    else
        V_bestRange_return_mps = V_bestRange_mps;
        Preq_return = Preq(idx_bestRange);
    end
    
    if strcmpi(assump.prop_type, 'electric')
        P_electric_outbound = Preq_outbound / (assump.eta_prop * assump.eta_electric);
        E_total_Wh = config.fuel_mass_kg * assump.battery_energy_density_Wh_per_kg;
        E_outbound_used_Wh = P_electric_outbound * t_outbound_hr;
        E_return_available_Wh = E_total_Wh - E_outbound_used_Wh;
        
        if E_return_available_Wh > 0 && P_electric_outbound > 0
            P_electric_return = Preq_return / (assump.eta_prop * assump.eta_electric);
            if P_electric_return > 0 && ~isnan(V_bestRange_return_mps)
                endurance_return_hr = E_return_available_Wh / P_electric_return;
                Range_return_km = V_bestRange_return_mps * endurance_return_hr * CONV.km_per_hr_per_mps;
                Range_return_nmi = Range_return_km * CONV.nmi_per_km;
                Range_mission_nmi = (Range_one_leg_req_km * CONV.nmi_per_km) + Range_return_nmi;
            end
        end
    else
        Pshaft_outbound_kW = (Preq_outbound / assump.eta_prop) / CONV.watts_per_kW;
        BSFC_kg_per_kWh = assump.BSFC_g_per_kWh / CONV.grams_per_kg;
        mdot_fuel_outbound_kg_per_hr = BSFC_kg_per_kWh * Pshaft_outbound_kW;
        fuel_outbound_used_kg = mdot_fuel_outbound_kg_per_hr * t_outbound_hr;
        fuel_return_available_kg = config.fuel_mass_kg - fuel_outbound_used_kg;
        
        if fuel_return_available_kg > 0 && mdot_fuel_outbound_kg_per_hr > 0
            Pshaft_return_kW = (Preq_return / assump.eta_prop) / CONV.watts_per_kW;
            mdot_fuel_return_kg_per_hr = BSFC_kg_per_kWh * Pshaft_return_kW;
            if mdot_fuel_return_kg_per_hr > 0 && ~isnan(V_bestRange_return_mps)
                endurance_return_hr = fuel_return_available_kg / mdot_fuel_return_kg_per_hr;
                Range_return_km = V_bestRange_return_mps * endurance_return_hr * CONV.km_per_hr_per_mps;
                Range_return_nmi = Range_return_km * CONV.nmi_per_km;
                Range_mission_nmi = (Range_one_leg_req_km * CONV.nmi_per_km) + Range_return_nmi;
            end
        end
    end
end

if isnan(Range_mission_nmi)
    if V_outbound_mps < min(V) || V_outbound_mps > max(V)
        warning('Mission range: 80 mph (%.2f m/s) outside velocity sweep range [%.2f, %.2f] m/s', ...
            V_outbound_mps, min(V), max(V));
    else
        warning('Mission range calculation failed - check energy/fuel calculations');
    end
end

[ROC_max_mps, idx_bestROC] = max(ROC);
ROC_max_fpm = ROC_max_mps * CONV.fpm_per_mps;
V_at_ROCmax_mph = V(idx_bestROC) * CONV.mph_per_mps;
gamma_deg_at_ROCmax = gamma_deg(idx_bestROC);

if V_req_top_mps >= min(V) && V_req_top_mps <= max(V)
    Preq_at_80mph = interp1(V, Preq, V_req_top_mps);
    ExcessP_at_80mph = interp1(V, ExcessP, V_req_top_mps);
else
    Preq_at_80mph = NaN;
    ExcessP_at_80mph = NaN;
end

Range_req_nmi_round_trip = requirements.Range_req_nmi_round_trip;
if ~isnan(Range_mission_nmi)
    if Range_mission_nmi >= Range_req_nmi_round_trip
        range_margin_pct = ((Range_mission_nmi / Range_req_nmi_round_trip) - 1) * 100;
        range_check = sprintf("PASS (+%.1f%%)", range_margin_pct);
    else
        range_margin_pct = (1 - (Range_mission_nmi / Range_req_nmi_round_trip)) * 100;
        range_check = sprintf("FAIL (-%.1f%%)", range_margin_pct);
    end
else
    range_check = "FAIL (cannot compute mission range)";
end

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


end
