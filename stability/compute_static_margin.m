function [SM, pass_fail] = compute_static_margin(x_cg_m, x_np_m, cbar_m, target_SM_min, target_SM_max)
% Computes static margin: SM = (x_np - x_cg) / cbar

SM = (x_np_m - x_cg_m) / cbar_m;
pass_fail = "PASS";
if SM < target_SM_min || SM > target_SM_max
    pass_fail = "FAIL";
end

end
