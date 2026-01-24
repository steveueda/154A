function post_convergence_summary(config, stability, PerfTable)
% POST CONVERGENCE SUMMARY
% Interprets and reports converged design results
% This function presents results only - no computations
%
% Inputs:
%   config - struct with converged geometry and weight:
%     .W - gross weight (N)
%     .Sref - wing reference area (m^2)
%     .cbar - mean aerodynamic chord (m)
%   stability - struct with stability analysis results:
%     .xCG - CG location (m from wing LE)
%     .xNP - neutral point location (m from wing LE)
%     .SM - static margin (fraction)
%   PerfTable - performance table from velocity_sweep.m

%% ====== CONFIGURATION SUMMARY ======
fprintf("\n========== FINAL CONFIGURATION SUMMARY ==========\n");
fprintf("Gross weight        : %.2f N (%.2f kg)\n", config.W, config.W/9.80665);
fprintf("Wing area           : %.3f m^2\n", config.Sref);
fprintf("Mean aero chord     : %.3f m\n", config.cbar);
fprintf("==================================================\n");

%% ====== STABILITY EXPLANATION ======
CG_percent = 100 * stability.xCG / config.cbar;
NP_percent = 100 * stability.xNP / config.cbar;
SM_percent = 100 * stability.SM;

fprintf("\n========== LONGITUDINAL STABILITY ==========\n");
fprintf("CG location         : %.3f m (%.1f %% MAC)\n", stability.xCG, CG_percent);
fprintf("Neutral point       : %.3f m (%.1f %% MAC)\n", stability.xNP, NP_percent);
fprintf("Static margin       : %.1f %% MAC\n", SM_percent);

% Stability verdict logic
if stability.SM > 0.05 && stability.SM < 0.15
    verdict = "POSITIVE static stability (within typical UAV range)";
elseif stability.SM <= 0
    verdict = "UNSTABLE (CG ahead of neutral point required)";
else
    verdict = "STABLE but overly conservative (excessive static margin)";
end

fprintf("Stability verdict   : %s\n", verdict);
fprintf("=============================================\n");

%% ====== VELOCITY SWEEP SUMMARY METRICS ======
% Extract key performance points
[~, idx_min_drag] = min(PerfTable.Drag_N);
[~, idx_max_ROC]  = max(PerfTable.ROC_mps);
[~, idx_max_rng]  = max(PerfTable.Range_km);

fprintf("\n========== PERFORMANCE HIGHLIGHTS ==========\n");
fprintf("Min drag speed      : %.2f m/s\n", PerfTable.V_mps(idx_min_drag));
fprintf("Max ROC speed       : %.2f m/s (ROC = %.2f m/s)\n", ...
        PerfTable.V_mps(idx_max_ROC), PerfTable.ROC_mps(idx_max_ROC));
fprintf("Max range speed     : %.2f m/s (Range = %.1f km)\n", ...
        PerfTable.V_mps(idx_max_rng), PerfTable.Range_km(idx_max_rng));
fprintf("=============================================\n");

%% ====== DISPLAY VELOCITY SWEEP TABLE ======
fprintf("\n========== VELOCITY SWEEP PERFORMANCE TABLE ==========\n");
disp(PerfTable);
fprintf("=======================================================\n\n");

end

