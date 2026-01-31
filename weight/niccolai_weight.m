function out = niccolai_weight(in)

W = in.W;
N = in.N;
S = in.S;
A = in.A;
Delta = in.Delta;
tr = in.tr;
tc = in.tc;
Ve = in.Ve;
lf = in.lf;
WF = in.WF;
D = in.D;
Sh = in.Sh;
bh = in.bh;
lh = in.lh;
hac = in.hac;
c = in.c;
ch = in.ch;
thr = in.thr;
Sv = in.Sv;
bv = in.bv;
tvr = in.tvr;
Weng = in.Weng;
Neng = in.Neng;

wing_weight_coeff = 96.948;
wing_weight_exp = 0.993;
wing_load_factor_exp = 0.65;
wing_AR_exp = 0.57;
wing_area_exp = 0.61;
wing_taper_exp = 0.36;
wing_velocity_exp = 0.5;
load_normalization = 1e5;
area_normalization = 100;
velocity_normalization_wing = 500;
velocity_normalization_fuse = 100;
length_normalization = 10;
taper_denominator = 2;

fuse_weight_coeff = 200;
fuse_weight_exp = 1.1;
fuse_load_factor_exp = 0.286;
fuse_length_exp = 0.857;
fuse_velocity_exp = 0.338;

ht_weight_coeff = 127;
ht_weight_exp = 0.458;
ht_load_factor_exp = 0.87;
ht_area_exp = 1.2;
ht_arm_exp = 0.483;
ht_thickness_exp = 0.5;

vt_weight_coeff = 98.5;
vt_weight_exp = 0.458;
vt_load_factor_exp = 0.87;
vt_area_exp = 1.2;
vt_thickness_exp = 0.5;
vt_factor = 2;
vt_area_factor = 0.5;
vt_span_factor = 0.5;

controls_weight_coeff = 1.066;
controls_weight_exp = 0.626;

prop_install_coeff = 2.575;
prop_install_exp = 0.922;

out.W_wing_lb = wing_weight_coeff * ((W*N/load_normalization)^wing_load_factor_exp * ...
    (A/cos(Delta))^wing_AR_exp * (S/area_normalization)^wing_area_exp * ...
    ((1+tr)/(taper_denominator*tc))^wing_taper_exp * (1+Ve/velocity_normalization_wing)^wing_velocity_exp)^wing_weight_exp;

out.W_fuse_lb = fuse_weight_coeff * ((W*N/load_normalization)^fuse_load_factor_exp * ...
    (lf/length_normalization)^fuse_length_exp * ((WF+D)/length_normalization) * ...
    (Ve/velocity_normalization_fuse)^fuse_velocity_exp)^fuse_weight_exp;

out.W_ht_lb = ht_weight_coeff * ((W*N/load_normalization)^ht_load_factor_exp * ...
    (Sh/area_normalization)^ht_area_exp * (lh/length_normalization)^ht_arm_exp * ...
    (bh/thr)^ht_thickness_exp)^ht_weight_exp;

out.W_vt_lb = vt_factor * vt_weight_coeff * ((W*N/load_normalization)^vt_load_factor_exp * ...
    (vt_area_factor*Sv/area_normalization)^vt_area_exp * ...
    (vt_span_factor*bv/tvr)^vt_thickness_exp)^vt_weight_exp;

out.W_controls_lb = controls_weight_coeff * W^controls_weight_exp;

out.W_prop_install_lb = prop_install_coeff * Weng^prop_install_exp * Neng;

out.W_struct_lb = out.W_wing_lb + out.W_fuse_lb + out.W_ht_lb + out.W_vt_lb;

end

