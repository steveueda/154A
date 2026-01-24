function nicIn = buildNiccolaiInputs(geomSI, W_guess_N, engine_mass_kg)
% Converts SI geometry and weight to US customary units for Niccolai weight estimation

conv = unit_conversions();

nicIn.W = W_guess_N * conv.lb_per_N;
nicIn.N = geomSI.n_ult;
nicIn.S = geomSI.Sref_m2 * conv.ft2_per_m2;
nicIn.A = geomSI.ARw;
nicIn.Delta = deg2rad(geomSI.sweep_w_deg);
nicIn.tr = geomSI.taper_w;
nicIn.tc = geomSI.tc_w;
nicIn.Ve = geomSI.Vmax_mps * conv.kt_per_mps;
nicIn.lf = geomSI.lf_m * conv.ft_per_m;
nicIn.WF = geomSI.wf_max_m * conv.ft_per_m;
nicIn.D = geomSI.df_max_m * conv.ft_per_m;
nicIn.Sh = geomSI.Sh_m2 * conv.ft2_per_m2;
nicIn.ARh = geomSI.ARh;
nicIn.bh = sqrt(geomSI.Sh_m2 * geomSI.ARh) * conv.ft_per_m;
nicIn.ch = sqrt(geomSI.Sh_m2 / geomSI.ARh) * conv.ft_per_m;
nicIn.thr = nicIn.ch * geomSI.tc_h * 12;
nicIn.c = sqrt(geomSI.Sref_m2 / geomSI.ARw) * conv.ft_per_m;
nicIn.lh = 0.5 * geomSI.lf_m * conv.ft_per_m;
nicIn.hac = 0.25;
nicIn.Sv = geomSI.Sv_m2 * conv.ft2_per_m2;
nicIn.ARv = geomSI.ARv;
nicIn.bv = sqrt(geomSI.Sv_m2 * geomSI.ARv) * conv.ft_per_m;
nicIn.cv = sqrt(geomSI.Sv_m2 / geomSI.ARv) * conv.ft_per_m;
nicIn.tvr = nicIn.cv * geomSI.tc_v * 12;
nicIn.Weng = engine_mass_kg * conv.lb_per_kg;
nicIn.Neng = 1;

end
