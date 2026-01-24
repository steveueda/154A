function x_np_m = compute_neutral_point(stab, Sref, ARw)
Vh = (stab.Sh_m2 * stab.lt_m) / (Sref * stab.cbar_m);
x_np_m = stab.x_ac_w_m + stab.cbar_m * (stab.eta_tail * (stab.a_t/stab.a_w) * (1 - stab.deps_dalpha) * Vh);
end
