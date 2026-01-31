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
geom.Sh_m2 = 0.035;         % Horizontal tail area (m^2)
geom.ARh = 4.0;            % Horizontal tail aspect ratio
geom.tc_h = 0.10;          % Horizontal tail thickness to chord ratio
geom.htail_xmc = 0.30;     % Horizontal tail max thickness location (fraction of chord)
geom.htail_sweep_deg = 0;  % Horizontal tail sweep (deg)
geom.Sv_m2 = 0.020;         % Vertical tail area (m^2)
geom.ARv = 1.8;            % Vertical tail aspect ratio
geom.tc_v = 0.10;          % Vertical tail thickness to chord ratio

%% Propulsion Assumptions
% Propulsion type: 'gas' or 'electric'
prop.type = 'electric';                 % 'gas' for gasoline engine, 'electric' for electric motor

% Common parameters
prop.eta_prop = 0.8;               % Propeller efficiency
prop.P_shaft_max_W = 375;          % Maximum shaft power (W)

% Gas engine parameters (used when prop.type = 'gas')
prop.BSFC_g_per_kWh = 500;         % Brake specific fuel consumption (g/kWh)
prop.fuel_density_kg_per_L = 0.74; % Fuel density (kg/L)

% Electric motor parameters (used when prop.type = 'electric')
prop.battery_energy_density_Wh_per_kg = 200;  % Battery energy density (Wh/kg) - typical LiPo: 150-250 Wh/kg
prop.battery_efficiency = 0.95;    % Battery discharge efficiency
prop.motor_controller_efficiency = 0.95;  % Motor controller efficiency
prop.eta_electric = prop.battery_efficiency * prop.motor_controller_efficiency;  % Overall electric efficiency

% Energy margins (applies to both gas and electric)
prop.reserve_frac = 0.1;          % Reserve energy fraction (10% reserve)
prop.tank_margin_frac = 0.10;      % Tank/battery sizing margin (10% margin)

%% Fixed Masses
% Set masses based on propulsion type
if strcmpi(prop.type, 'electric')
    mass.engine_mass_kg = 0.31;  % Electric motor + controller (Motor is 0.235 kg, controller is approx 0.075 kg)
    mass.fuel_system_mass_kg = 0.15 * 0.40;  % Battery management system, wiring, cooling 
else
    mass.engine_mass_kg = 2.1;  % Gas engine
    mass.fuel_system_mass_kg = 0.20 * 0.40;  % Fuel system (pump, lines, etc.) 
end

mass.avionics_mass_kg = 0.3;  
mass.landing_gear_mass_kg = 0.0; % I set at zero for now as Josh suggested we dont need one

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
% x = 0 at front tip of plane (m)
cg.wing_LE_from_front_m = 0.40;  % Distance from front tip to wing leading edge (m)
cg.payload_x_m = 0.535;        % Payload CG location (m) at the cg pretty much so that cg doesnt change when dropped
cg.engine_x_m = 0.1;         % Engine/motor CG location (m) at the front of the plane 
cg.battery_x_m = 0.48;        % Battery CG location (m) near the cg beccause heavy
cg.avionics_x_m = 1.1;       % Avionics CG location (m) towards the back
cg.landing_gear_x_m = 0;   % Landing gear CG location (m) (Not used)
cg.fuel_x_m = 0.62;           % Gas fuel tank CG location (m) (not used for electric)
cg.fuel_system_x_m = 0.60;    % Fuel system/BMS CG location (m) (not used for electric)

%% Weight Correction Factors
% Correction factors for Niccolai weight estimates (applied to calculated weights only)
% User input values (avionics, fuel_system, payload, battery) are already corrected in their definitions
k_corr.wing        = 2.00;
k_corr.fuselage    = 1.25;
k_corr.htail       = 0.70;
k_corr.vtail       = 0.60;
k_corr.controls    = 0.11;
k_corr.prop_install = 0.55;

%% Analysis Controls
% Iteration settings
ctrl.max_iter = 8;
ctrl.tol_frac = 0.01;      % 1% relative convergence tolerance

% Velocity sweep for performance analysis
ctrl.V_min_fraction = 0.4;  % Minimum velocity as fraction of stall speed (e.g., 0.4 = 40% of stall)
ctrl.V_max_fraction = 1.15;  % Maximum velocity as fraction of cruise speed (e.g., 2.0 = 200% of cruise)

% Plot x-axis limits (as fraction of velocity range)
ctrl.power_plot_xlim_low = 0;   % Power plot lower limit (0.85 = 15% below min velocity)
ctrl.power_plot_xlim_high = 1.1;   % Power plot upper limit (1.1 = 10% beyond max velocity)
ctrl.drag_plot_xlim_low = 1;    % Drag plot lower limit (0.85 = 15% below min velocity)
ctrl.drag_plot_xlim_high = 1.1;    % Drag plot upper limit (1.1 = 10% beyond max velocity)

% Velocity label positions (NaN = use relative offset from marker, otherwise use absolute coordinates)
% Power vs. Velocity plot labels
ctrl.label_power.v_stall_x = 17;      % x-coordinate (m/s) for v_stall label (NaN = relative offset)
ctrl.label_power.v_stall_y = NaN;      % y-coordinate (W) for v_stall label (NaN = relative offset)
ctrl.label_power.v_minpower_x = 12;   % x-coordinate (m/s) for v_minpower label (NaN = relative offset)
ctrl.label_power.v_minpower_y = 76;    % y-coordinate (W) for v_minpower label
ctrl.label_power.v_range_x = 20.3;       % x-coordinate (m/s) for v_range label
ctrl.label_power.v_range_y = NaN;      % y-coordinate (W) for v_range label (NaN = relative offset)
ctrl.label_power.v_max_x = NaN;        % x-coordinate (m/s) for v_max label (NaN = relative offset)
ctrl.label_power.v_max_y = NaN;        % y-coordinate (W) for v_max label (NaN = relative offset)

% Drag vs. Velocity plot labels
ctrl.label_drag.v_mindrag_x = 20.2;     % x-coordinate (m/s) for v_mindrag label (NaN = relative offset)
ctrl.label_drag.v_mindrag_y = 4.7;     % y-coordinate (N) for v_mindrag label (NaN = relative offset)
ctrl.label_drag.v_stall_x = 18;       % x-coordinate (m/s) for v_stall label (NaN = relative offset)
ctrl.label_drag.v_stall_y = 4.2;       % y-coordinate (N) for v_stall label (NaN = relative offset)


