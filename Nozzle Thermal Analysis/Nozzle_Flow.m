function [h, T_aw, i] = Nozzle_Flow(machPosition, machDiameter, convergingAngle, divergingAngle, curveRadius, ...
    stagnationTemperature, stagnationPressure, stagnationDensity, CEAgamma, CEAcp, position, fzeroMargin)
% Calculates heat flux coefficient (h) along a nozzle, temperature of the
% gas along the nozzle wall (T_aw), and the index in the position vector
% corresponding to the throat (i)
% All inputs must be in mm, see below for the default values
% Uses pipe-flow equation for fully developed flow (see pdf in code folder)
if nargin==0
    machPosition = 27.72; % mm
    machDiameter = 14.56; % mm
    convergingAngle = pi/4; % angle of converging section in rad
    divergingAngle = pi/12; % angle of diverging section in rad
    curveRadius = 0; % mm
    stagnationTemperature = 2746; % K
    stagnationPressure = 40e5; % Pa
    stagnationDensity = 4.754; % kg/m^3
    CEAgamma = [[0;27.72;64.12], [1.2071;1.2259;1.2556]]; % [position in mm], [gamma]
    CEAcp = [[0;27.72;64.12], [1.9263e3;1.7067e3;1.5e3]]; % [position in mm], [cp]
    position = [0, 5, 10, 15, 20, 25, 26, 27, 27.5 27.72, 28, 29, 30, 35, 40, 45, 50, 55, 60, 64.12] % positions in mm
    fzeroMargin = 100; % margin of boundary used in MachNumber
end

close all
gamma = Gamma(CEAgamma, position);
cp = Cp(CEAcp, position);
% diameter = Diameter(convergingAngle, divergingAngle, machPosition, machDiameter, position);
diameter = DiameterCurved(convergingAngle, divergingAngle, curveRadius, machPosition, machDiameter, position);
areaRatio = AreaRatio(diameter, machDiameter);
supersonic = position>machPosition;
machNumber = MachNumber(areaRatio, supersonic, gamma, fzeroMargin);
tempRatio = TempRatio(machNumber, gamma);
temp = tempRatio .* stagnationTemperature;
mu = DynamicViscosity(temp);
presRatio = PresRatio(machNumber, gamma);
pres = presRatio .* stagnationPressure;
rhoRatio = RhoRatio(machNumber, gamma);
rho = rhoRatio .* stagnationDensity;
Pr = Prandtl(temp);
u = Velocity(machNumber, gamma, pres, rho);
diameter_m = diameter ./ 1000; % converting from mm to m

Re = rho .* u .* diameter_m ./ mu
St = 0.023 .* Re.^(-0.2) ./ Pr.^0.6;
h = St .* rho .* u .* cp %should be in (J/(s*m^2*K))
% massFlux = rho .* u
figure
plot(position, diameter/2)
ylabel("radius (mm)")
xlabel("position (mm)")
title("Nozzle shape")

figure
plot(position, h)
xlabel("Position (mm)")
ylabel("heat flux coefficient h (J/(s*m^2*K))")
title("Heat flux coefficient along the length of the nozzle")

T_aw_Ratio = T_awRatio(Pr, gamma, machNumber);
T_aw = T_aw_Ratio .* stagnationTemperature
i = find(position==machPosition);

figure
plot(position, T_aw)
xlabel("Position (mm)")
ylabel("adiabatic wall temperature (K)")
title("Adiabatic wall temperature along the length of the nozzle")

% Semi-infinite wall heat flow analysis
% T_wall0 = 298;
% c_p_graphite = 1.534e6; % J/(m^3*K)
% d = 0.01; % m, 0.01 m = 10 mm
% C_1 = -c_p_graphite*d/h(i)*log(T_aw(i)-T_wall0)
% time = [0,5,10,15,20,25,30];
% T_wall = T_aw(i)-exp(-h(i)./(c_p_graphite*d).*(time+C_1));
% figure
% plot(time, T_wall)
% xlabel("time (s)")
% ylabel("10 mm wall temperature (K)")