function Pr = Prandtl (temp)
% Experimental data for Air at 10 atm (see pdf in code folder)
% Linear regression from (1300 K, .705) to (2000 K, .7)
% Quadratic regression from (2000 K, .7) to (2500 K, .684) and to 2750 K
% Does not depend on nozzle design, so no user input required
Pr_lin = -6.786e-6 .* temp + 0.7138;
Pr_quad = -4.107e-8 .* temp.^2 + 1.534e-4 .* temp + 0.5573;
Pr = Pr_lin .* (temp<=2000) + Pr_quad .* (temp>2000);