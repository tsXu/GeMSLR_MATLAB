function [ Q, H, m2, tits ] = subspace_iteration_gemslr(PRE, levi, k, its, orthtol, resits)
%[ Q, H, m2, tits ] = subspace_iteration_gemslr(PRE, levi, k, its, orthtol, resits)
%   Use subspace iteration to generate Ritz values and vectors

n = size(PRE.Levs{levi}.C,1);

% initial guess
B = rand(n,k);

% main loop
tits = 0;
for i = 1:its
    tits = tits + k;
    [Q, ~] = qr(B,0);
    B = EBinvFCinv(PRE, Q, levi, resits);
    %Ynorm = norm(B-Q*Q'*B,'fro');
    %if(Ynorm < orthtol)
    %    break;
    %end
end

%[Q, ~] = qr(B,0);
H = Q'*B;

%[Q1, H] = schur(H);
%Q = Q*Q1;

m2 = k;

end