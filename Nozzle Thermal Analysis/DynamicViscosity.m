function mu = DynamicViscosity (temp)
% Calculates mu (Dynamic Viscosity) using Sutherland's Law for Air
% (see pdf in code folder)
% Does not depend on nozzle design, so no user input required
mu_0 = 1.716e-5;
temp_0 = 273.15;
S = 110.4;
mu = mu_0.*(temp./temp_0).^(3/2).*(temp_0+S)./(temp+S);