function Pres_Ratio = PresRatio (Mach_Number, gamma)
% Calculates pressure ratio from the mach number and 
% gamma (specific heat ratio) of the flow. Uses the
% isentropic flow relation.
exponent = -gamma ./ (gamma -1);
Pres_Ratio = (1+(gamma-1)./2 .* Mach_Number.^2).^ exponent;