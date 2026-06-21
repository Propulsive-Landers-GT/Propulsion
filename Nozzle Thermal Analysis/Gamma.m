function gamma = Gamma (CEAgamma, position)
% calculate gamma (heat capacity ratio) at chosen positions by linearly
% interpolating data from CEA. If a single value is provided, gamma will
% be calculated as constant
% for the original nozzle design, gamma = 1.206 + 7.5943e-4 .* position;
linReg = [ones(length(CEAgamma(:,1)),1), CEAgamma(:,1)]\CEAgamma(:,2);
intercept = linReg(1);
slope = linReg(2);
gamma = intercept + slope .* position;