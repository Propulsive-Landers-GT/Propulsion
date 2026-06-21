function cp = Cp (CEAcp, position)
% calculate cp (specific heat) at chosen positions by linearly
% interpolating data from CEA. If a single value is provided, cp will
% be calculated as constant
% for the original nozzle design, cp = 1913 - 6.599 .* position;;
linReg = [ones(length(CEAcp(:,1)),1), CEAcp(:,1)]\CEAcp(:,2);
intercept = linReg(1);
slope = linReg(2);
cp = intercept + slope .* position;
