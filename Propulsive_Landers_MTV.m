% Propulsive Landers MTV - Valve Flow Analysis with user-specified ΔP
% Units: Final flow rates in kg/s

% Assumptions made in this script [READ THIS PLEASE]
% - Cv values provided are valid and based on manufacturer data for water at 60°F.
% - Cv is assumed to be a property of the valve only (geometry-based), independent of the working fluid.
% - Liquid nitrous oxide (N2O) is treated as a single-phase incompressible liquid with constant density = 745 kg/m^3 at ~20 °C.
% - Specific gravity (SG) is computed relative to water at 1000 kg/m^3.
% - Pressure drop across the valve (ΔP) is supplied by the user as a constant value.
% - ΔP is assumed to remain above N2O vapor pressure at operating temperature 
%   (no flashing, cavitation, or two-phase effects modeled).
% - Cv equation used: Q [GPM] = Cv * sqrt(ΔP [psi] / SG).
% - Conversion from GPM → m^3/s → kg/s assumes constant fluid density (no thermal effects).
% - No line losses, injector drops, or transient dynamics are included in this model.
% - This is a steady-state, sizing-level analysis for valve selection only.

% Inputs 
rho = 745;             % Density of liquid N2O [kg/m^3] (example ~20 °C)
SG  = rho / 1000;      % Specific gravity relative to water

DeltaP_bar = 2;        % <-- YOU enter ΔP across the valve [bar]
DeltaP_psi = DeltaP_bar * 14.5038;   % Convert to psi for Cv equation

% Conversion constants
GPM_to_m3s = 6.309e-5; % 1 GPM = 6.309e-5 m^3/s

% Valve Data
ValveOpen = [0 10 20 30 40 50 60 70 80 90 100];

CV30 = [0.000 0.000 0.100 0.172 0.324 0.429 0.649 0.873 1.350 1.749 2.435];
CV60 = [0.000 0.000 0.120 0.236 0.539 0.643 1.081 1.587 2.615 3.664 5.525];
CV90 = [0.000 0.100 0.200 0.400 0.600 0.800 1.500 2.200 3.800 5.400 6.900];

% Flow Rate Calculation 
% Cv equation: Q_gpm = Cv * sqrt(ΔP / SG)
Q30_gpm = CV30 .* sqrt(DeltaP_psi ./ SG);
Q60_gpm = CV60 .* sqrt(DeltaP_psi ./ SG);
Q90_gpm = CV90 .* sqrt(DeltaP_psi ./ SG);

% Convert GPM → m^3/s → kg/s
Q30 = Q30_gpm * GPM_to_m3s * rho;
Q60 = Q60_gpm * GPM_to_m3s * rho;
Q90 = Q90_gpm * GPM_to_m3s * rho;

figure;

% Subplot 1: 30 deg
subplot(2,2,1);
plot(ValveOpen, Q30, '-o', 'LineWidth', 1.5);
xlabel('Valve Opening (%)'); ylabel('Flow Rate (kg/s)');
title('30° Valve'); grid on;

% Subplot 2: 60 deg
subplot(2,2,2);
plot(ValveOpen, Q60, '-o', 'LineWidth', 1.5);
xlabel('Valve Opening (%)'); ylabel('Flow Rate (kg/s)');
title('60° Valve'); grid on;

% Subplot 3: 90 deg
subplot(2,2,3);
plot(ValveOpen, Q90, '-o', 'LineWidth', 1.5);
xlabel('Valve Opening (%)'); ylabel('Flow Rate (kg/s)');
title('90° Valve'); grid on;

% Subplot 4: Combined
subplot(2,2,4);
plot(ValveOpen, Q30, '-o', 'LineWidth', 1.5); hold on;
plot(ValveOpen, Q60, '-s', 'LineWidth', 1.5);
plot(ValveOpen, Q90, '-^', 'LineWidth', 1.5);
xlabel('Valve Opening (%)'); ylabel('Flow Rate (kg/s)');
title('Comparison of All Valves');
legend('30°','60°','90°','Location','northwest'); grid on;

fprintf('Valve ΔP used = %.2f bar (%.2f psi)\n', DeltaP_bar, DeltaP_psi);
