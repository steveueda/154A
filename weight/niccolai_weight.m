function out = niccolai_weight(in)
% NICCOLAI WEIGHT ESTIMATION
% Computes component weights using Niccolai formulas (exact equations)
% 
% Input: struct 'in' with fields (all in US customary units)
% Output: struct 'out' with component weights (lb)

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
cv = in.cv;
tvr = in.tvr;
Weng = in.Weng;
Neng = in.Neng;

% Wing Weight (Niccolai) - EXACT EQUATION
out.W_wing_lb = 96.948*((W*N/10^5)^0.65*(A/cos(Delta))^0.57*(S/100)^0.61*((1+tr)/(2*tc))^0.36*(1+Ve/500)^0.5)^0.993;

% Fuselage Weight (Niccolai) - EXACT EQUATION
out.W_fuse_lb = 200*((W*N/10^5)^0.286*(lf/10)^0.857*((WF+D)/10)*(Ve/100)^0.338)^1.1;

% Horizontal Tail Weight (Niccolai) - EXACT EQUATION
out.W_ht_lb = 127*((W*N/10^5)^0.87*(Sh/100)^1.2*(lh/10)^0.483*(bh/thr)^0.5)^0.458;

% Vertical Tail Weight (Niccolai) - EXACT EQUATION
out.W_vt_lb = (2)*98.5*((W*N/10^5)^0.87*((.5)*Sv/100)^1.2*((.5)*bv/tvr)^0.5)^0.458;

% Surface Controls Weight (Niccolai) - EXACT EQUATION
out.W_controls_lb = 1.066*W^0.626;

% Propulsion Installation Weight (Niccolai) - EXACT EQUATION
out.W_prop_install_lb = 2.575*(Weng)^0.922*Neng;

% Total Structural Weight
out.W_struct_lb = out.W_wing_lb + out.W_fuse_lb + out.W_ht_lb + out.W_vt_lb;

end

