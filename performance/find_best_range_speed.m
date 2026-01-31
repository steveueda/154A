function [V_bestRange_mps, Preq_bestRange_W] = find_best_range_speed(W_N, Sref_m2, components, env, aero, dragParams, assump, V_min_mps, V_max_mps)
% Finds the best range speed (maximizes V/Preq) for a given weight and geometry
%
% Inputs:
%   W_N - Aircraft weight (N)
%   Sref_m2 - Wing reference area (m^2)
%   components - Component struct array for drag calculation
%   env - Environment struct (rho, mu, a)
%   aero - Aerodynamics struct (ARw, ew)
%   dragParams - Drag parameters struct
%   assump - Assumptions struct
%   V_min_mps - Minimum velocity to search (m/s)
%   V_max_mps - Maximum velocity to search (m/s)
%
% Outputs:
%   V_bestRange_mps - Best range speed (m/s)
%   Preq_bestRange_W - Power required at best range speed (W)

N_sweep_points = 50;
fallback_speed_fraction = 0.7;

if V_min_mps >= V_max_mps
    V_bestRange_mps = V_max_mps * fallback_speed_fraction;
    drag = compute_drag(V_bestRange_mps, W_N, Sref_m2, components, env, aero, dragParams);
    power = compute_power(drag, V_bestRange_mps, assump);
    Preq_bestRange_W = power.Preq_W;
    return;
end

V = linspace(V_min_mps, V_max_mps, N_sweep_points)';
Preq = zeros(size(V));

for i = 1:length(V)
    drag = compute_drag(V(i), W_N, Sref_m2, components, env, aero, dragParams);
    power = compute_power(drag, V(i), assump);
    Preq(i) = power.Preq_W;
end

valid_indices = (Preq > 0);
if any(valid_indices)
    range_metric = V(valid_indices) ./ Preq(valid_indices);
    [~, idx_local_max] = max(range_metric);
    V_bestRange_mps = V(valid_indices(idx_local_max));
    Preq_bestRange_W = Preq(valid_indices(idx_local_max));
else
    V_bestRange_mps = NaN;
    Preq_bestRange_W = NaN;
    warning('Could not find a valid best range speed (Preq was not positive).');
end

end

