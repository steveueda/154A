function power = compute_power(drag, V, assump)
% Computes power required and available

power.Preq_W = drag.D_N * V;
power.Pavail_W = assump.eta_prop * assump.P_shaft_max_W;
power.ExcessP_W = power.Pavail_W - power.Preq_W;

end
