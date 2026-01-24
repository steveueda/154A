function [x_cg_m, total_mass_kg, CGTable] = compute_cg(compMassLoc)
% Computes CG location from component masses and positions

n = length(compMassLoc);
names = strings(n,1);
masses = zeros(n,1);
x_pos = zeros(n,1);

for i = 1:n
    names(i) = string(compMassLoc(i).name);
    masses(i) = compMassLoc(i).mass_kg;
    x_pos(i) = compMassLoc(i).x_m;
    if masses(i) < 0
        error("Component '%s' has negative mass: %.3f kg", names(i), masses(i));
    end
end

total_mass_kg = sum(masses);
if total_mass_kg <= 0
    error("Total mass is zero or negative: %.3f kg", total_mass_kg);
end

moments = masses .* x_pos;
x_cg_m = sum(moments) / total_mass_kg;
CGTable = table(names, masses, x_pos, moments, ...
    'VariableNames', {'name', 'mass_kg', 'x_m', 'moment_kgm'});

end
