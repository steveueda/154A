function Sref_m2 = size_wing(W_N, rho, V_stall, CL_max)
Sref_m2 = W_N / (0.5 * rho * V_stall^2 * CL_max);
end
