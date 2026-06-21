function Diameter = Diameter(convergingAngle, divergingAngle, machPosition, machDiameter, position)
% Calculates diameter (in mm) along the nozzle. Converging section goes
% from position = 0 to machPosition. Diverging section goes from 
% position = machPosition till final position value. Does not
% handle curved nozzle throats.
diameter1 = (machDiameter - 2*tan(convergingAngle).*(position - machPosition)).*(position <= machPosition); % Converging diameter
diameter2 = (machDiameter + 2*tan(divergingAngle).*(position - machPosition)).*(position > machPosition); % Diverging diameter
Diameter = diameter1 + diameter2;
% plot(position, Diameter)