function [ L, W, D ] = unweighted_laplacian( A )
%[ L, W, D ] = unweighted_laplacian( A )
%   Create unweighted graph Laplacian of |A|+|A'|

n = size(A,1);

% make symmetric and remove diagonal
A = abs(A);
A = A + A';

% remove diagonal
A(1:n+1:end) = 0;

% generate laplacian
[ii, jj, ~] = find(A);
nnz= length(ii);
aa = ones(nnz, 1);

W = sparse(ii, jj, aa, n, n);
D = sum(W,2);
D = spdiags(D,0,n,n);
L = D - W;

end

