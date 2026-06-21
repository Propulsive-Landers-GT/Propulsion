function Temp_Ratio = TempRatio (Mach_Number, gamma)
% Calculates temperature ratio from the mach number and 
% gamma (specific heat ratio) of the flow. Uses the
% isentropic flow relation.
Temp_Ratio = 1 ./ (1+(gamma-1)./2 .* Mach_Number.^2);