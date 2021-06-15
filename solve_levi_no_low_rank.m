function [ y ] = solve_levi_no_low_rank( PRE, x, levi )
%solve_levi_no_low_rank( PRE, x, levi )
%   Apply a solve with levi matrix without low-rank correction
%   For example, solve with A is when levi == 1
%   Solve with C_i at level i is when levi == i+1 (C_i is the matrix on level levi+1)
%   PRE:    Preconditioning data structure
%   x:      The input vector
%   levi:   Level to apply the solve

nlev = PRE.nlev;

if levi == nlev
    % returen on the last level
    y = ( PRE.Levs{levi}.UB \ ( PRE.Levs{levi}.LB \ x ) );
    return;
end

nB = size(PRE.Levs{levi}.LB, 1);
% split
xB = x(1:nB,:);
xC = x(nB+1:end,:);
% solve with B
zB = PRE.Levs{levi}.UB \ ( PRE.Levs{levi}.LB \ xB );
zC = xC - PRE.Levs{levi}.E * zB;
% solve with C
yC = solve_levi_no_low_rank( PRE, zC, levi + 1 );
%yC = 1/(1-PRE.theta)*(PRE.UC\(PRE.LC\zC)) + PRE.Z*(PRE.G*(PRE.Y'*zC));
% back to last level
yB = zB - PRE.Levs{levi}.UB \ ( PRE.Levs{levi}.LB \ (PRE.Levs{levi}.F * yC));
y = [yB; yC];

end

