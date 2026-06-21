function Diameter = DiameterCurved(convergingAngle, divergingAngle, curveRadius, machPosition, machDiameter, position)
% Calculates diameter (in mm) along the nozzle. Converging section goes
% from position = 0 to leftTangentPosition. Curved section goes from 
% leftTangentPosition to rightTangentPosition. Diverging section goes 
% from rightTangentPosition till final position value.
leftTangentPosition = machPosition - curveRadius*sin(convergingAngle);
rightTangentPosition = machPosition + curveRadius*sin(divergingAngle);
leftTangentDiameter = machDiameter + 2*(curveRadius-(curveRadius^2-(machPosition-leftTangentPosition)^2)^.5);
rightTangentDiameter = machDiameter + 2*(curveRadius-(curveRadius^2-(machPosition-rightTangentPosition)^2)^.5);
diameter1 = (leftTangentDiameter - 2*tan(convergingAngle).*(position - leftTangentPosition)).*(position <= leftTangentPosition); % Converging diameter
diameter2 = (machDiameter + 2.*(curveRadius-(curveRadius^2-(machPosition-position).^2).^.5)).*(position > leftTangentPosition).*(position <= rightTangentPosition); % curve diameter
diameter3 = (rightTangentDiameter + 2*tan(divergingAngle).*(position - rightTangentPosition)).*(position > rightTangentPosition); % Diverging diameter
Diameter = diameter1 + diameter2 + diameter3;
