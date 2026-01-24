%% USER INPUTS

%% Mission Requirements
req.Range_req_km = 30 * 1852 / 1000;  % Round trip range (km) 
req.V_stall_mps = 40 * 0.44704;       % Stall speed requirement (m/s) 
req.V_design_mps = 80 * 0.44704;       % Design cruise/top speed (m/s) 
req.payload_mass_kg = 4.5 * 0.45359237;  % Payload mass (kg) - 4.5 lb

%% Aerodynamic Assumptions
aero.CL_max = 1.5;                 % Maximum lift coefficient (for stall sizing)
aero.ew = 0.85;                    % Wing Oswald efficiency factor

aero.Q = 1.05;                     % Interference factor 
aero.CD_misc = 0.0005;             % Miscellaneous drag coefficient
aero.CD_LP = 0.0010;               % Leakage and protuberance drag coefficient

%% Geometry Assumptions
% Flight conditions 
geom.rho = 1.225;          % Air density (kg/m^3)
geom.mu  = 1.789e-5;       % Dynamic viscosity (Pa s)
geom.a   = 340;            % Speed of sound (m/s)

% Fuselage dimensions
geom.lf_m = 1.50;          % Fuselage length (m)
geom.wf_max_m = 0.20;      % Fuselage max width (m)
geom.df_max_m = 0.20;      % Fuselage max depth (m)

% Wing parameters
geom.ARw = 7.9;            % Wing aspect ratio
geom.taper_w = 0.60;       % Wing taper ratio
geom.tc_w = 0.12;          % Wing thickness to chord ratio
geom.sweep_w_deg = 0;      % Wing sweep at quarter chord (deg)
geom.wing_xmc = 0.30;      % Wing max thickness location (fraction of chord)

% Tail parameters
geom.Sh_m2 = 0.10;         % Horizontal tail area (m^2)
geom.ARh = 4.0;            % Horizontal tail aspect ratio
geom.tc_h = 0.10;          % Horizontal tail thickness to chord ratio
geom.htail_xmc = 0.30;     % Horizontal tail max thickness location (fraction of chord)
geom.htail_sweep_deg = 0;  % Horizontal tail sweep (deg)
geom.Sv_m2 = 0.08;         % Vertical tail area (m^2)
geom.ARv = 1.8;            % Vertical tail aspect ratio
geom.tc_v = 0.10;          % Vertical tail thickness to chord ratio


%% Propulsion Assumptions
prop.eta_prop = 0.8;               % Propeller efficiency
prop.P_shaft_max_W = 1000;         % Maximum shaft power (W)
prop.BSFC_g_per_kWh = 500;         % Brake specific fuel consumption (g/kWh)
prop.fuel_density_kg_per_L = 0.74; % Fuel density (kg/L)

% Fuel margins
prop.reserve_frac = 0.15;          % Reserve fuel fraction
prop.tank_margin_frac = 0.10;      % Tank sizing margin

%% Fixed Masses
mass.engine_mass_kg = 2.1;
mass.avionics_mass_kg = 0.45;
mass.landing_gear_mass_kg = 0.70;
mass.fuel_system_mass_kg = 0.20;

mass.W0_N = 216.7;         % Initial gross weight guess (N)
mass.n_ult = 5.7;          % Ultimate load factor

%% Stability Parameters
stab.eta_tail = 0.95;      % Tail dynamic pressure ratio
stab.a_w = 2*pi;           % Wing lift curve slope (1/rad)
stab.a_t = 2*pi;           % Tail lift curve slope (1/rad)
stab.deps_dalpha = 0.30;   % Downwash gradient
stab.lt_frac = 0.5;        % Tail arm as fraction of fuselage length

% Static margin targets
stab.SM_min = 0.08;        % Minimum acceptable static margin 
stab.SM_max = 0.15;        % Maximum acceptable static margin 

%% Component Locations for CG Analysis
% x = 0 at wing leading edge (m)
cg.payload_x_m = 0.30;
cg.engine_x_m = 0.10;
cg.avionics_x_m = 0.25;
cg.landing_gear_x_m = 0.15;
cg.fuel_x_m = 0.22;        % Fuel CG location (m from wing LE)
cg.fuel_system_x_m = 0.20;

%% Analysis Controls
% Iteration settings
ctrl.max_iter = 8;
ctrl.tol_frac = 0.01;      % 1% relative convergence tolerance

% Velocity sweep for performance analysis
ctrl.Vmin = 10;            % m/s
ctrl.Vmax = 45;            % m/s
ctrl.N_sweep = 30;         % Number of velocity points


