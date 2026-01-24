function drag = compute_drag(V, W, Sref, components, env, aero, dragParams)
% Computes parasite and induced drag at a single velocity

M = V / env.a;
CDp_sum = 0;

for i = 1:numel(components)
    c = components(i);
    Re = env.rho * V * c.Lref / env.mu;
    Cf = 0.455 / ((log10(Re)^2.58) * ((1 + 0.144*M^2)^0.65));
    
    if strcmpi(c.type,"wing")
        K = (1 + (0.6/c.xmc)*c.tc + 100*(c.tc^4)) * (1.34*(M^0.18)*(cos(c.sweep_m_rad)^0.28));
    elseif strcmpi(c.type,"fuselage")
        K = 1 + 60/(c.fineness^3) + c.fineness/400;
    else
        error("Unknown component type: %s", c.type);
    end
    
    CDp_sum = CDp_sum + Cf * K * dragParams.Q * (c.Swet / Sref);
end

drag.CDp = CDp_sum + dragParams.CD_misc + dragParams.CD_LP;
CL = 2*W/(env.rho*V^2*Sref);
drag.CDi = CL^2 / (pi * aero.ARw * aero.ew);
drag.CDtotal = drag.CDp + drag.CDi;
drag.D_N = 0.5*env.rho*V^2*Sref*drag.CDtotal;

end
