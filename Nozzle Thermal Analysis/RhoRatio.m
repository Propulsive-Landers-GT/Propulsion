function Rho_Ratio = RhoRatio (Mach_Number, gamma)
% Calculates rho (density) ratio from the mach number and 
% gamma (specific heat ratio) of the flow. Uses the
% isentropic flow relation.
exponent = -1 ./ (gamma -1);
Rho_Ratio = (1+(gamma-1)./2 .* Mach_Number.^2).^ exponent;