function T_aw_Ratio = T_awRatio(Pr, gamma, machNumber)
% Calculates the ratio between gas temperature at the nozzle wall and the
% stagnation temperature from the Prandtl number, gamma (specific heat capacity),
% and mach number (see pdf in code folder)
gamma_frac = (gamma - 1)./2;
T_aw_Ratio = (1+Pr.^(1/3).*gamma_frac.*machNumber.^2)./(1+gamma_frac.*machNumber.^2);