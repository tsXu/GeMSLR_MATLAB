function y = twoLevSol(PRE, x)
nB = size(PRE.LB);
%
xB = x(1:nB);
xC = x(nB+1:end);
%
zB = PRE.UB \ ( PRE.LB \ xB );
zC = xC - PRE.ET * zB;
%
zC = zC + PRE.Z*(PRE.G*(PRE.Y'*zC));
yC = (PRE.UC\(PRE.LC\zC));
%yC = 1/(1-PRE.theta)*(PRE.UC\(PRE.LC\zC)) + PRE.Z*(PRE.G*(PRE.Y'*zC));
%
yB = zB - PRE.UB \ ( PRE.LB \ (PRE.F * yC));
y = [yB; yC];
end