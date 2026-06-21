function Mach_Number = MachNumber(Area_Ratio, Supersonic, gamma, margin)
% Calculates mach number from area ratio, a bool flag whether the flow is
% supersonic, and the gamma (specific heat ratio) of the flow. Uses the
% isentropic flow relation. Since the relation can't be solved explicitly,
% this uses the fzero function. If the numbers produced are inconsistent,
% try increasing margin from the default of 100.
Mach_Number = [];
for i = [Area_Ratio; Supersonic; gamma]
    % i
    if abs(i(1)-1) < 1e-3
        Mach_Number = [Mach_Number, 1];
    else
        gp1 = i(3) + 1;
        gm1 = i(3) - 1;
        C1 = (gp1/2).^(-(gp1)./(2*gm1));
        C2 = gm1/2;
        C3 = gp1./2./gm1;
        AreaRatioZero = @(M) C1.*(1+C2.*M.^2).^C3/M - i(1);
        Mach_Number = [Mach_Number, fzero(AreaRatioZero,[1/margin,1]*margin.^i(2))];
    end
end