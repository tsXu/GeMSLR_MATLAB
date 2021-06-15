function [ y ] = SCinv( PRE, x, levi )
%[ y ] = SCinv( PRE, x, levi )
%   Approximately solve I - SC^{-1};
%   PRE:    Preconditioning data structure
%   x:      The input vector
%   levi:   Level to apply the solve

y = x;

% compute Cinv*x
y1 = solve_levi(PRE, x, levi+1);

% compute (C - EB^{-1}F) y1
y = y - PRE.Levs{levi}.C * y1;

% compute Fy1
y1 = PRE.Levs{levi}.F*y1;
% compute BinvFy1
y1 = ( PRE.Levs{levi}.UB \ ( PRE.Levs{levi}.LB \ y1 ) );

% compute EBinvFy1
y = y + PRE.Levs{levi}.E*y1;

end

