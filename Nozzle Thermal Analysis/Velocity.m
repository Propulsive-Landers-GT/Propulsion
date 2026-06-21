function u = Velocity (mach, gamma, pres, dens)
% calculates gas flow velocity by finding the speed of sound and
% multiplying by the mach number
a = (gamma .* pres ./ dens).^0.5;
u = a .* mach;