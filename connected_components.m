function [ comps ] = connected_components( A )
%[ comps ] = connected_components( A )
%   find connected comoonents in A
%   A: matrix
%   comps: cell, comps{i} is the ith connected component

n = size(A,1);

ncomp = 0;
marker = zeros(1,n); % helper array

for i = 1:n
    if marker(i) == 0
        % unvisited, visit from here
        [marker, compi] = visit(A, i, marker);
        ncomp = ncomp + 1;
        comps{ncomp} = compi;
    end
end

end

function [marker, compi] = visit(A, i, marker)
    % visit from node i
    n = size(A,1);
    q = zeros(1,n);
    qs = 1;
    qe = 1;
    q(1) = i;
    marker(i) = 1;
    
    while qs <= qe
        idx = q(qs);
        [ii,~,~] = find(A(:,idx));
        ii = ii(marker(ii)==0);
        ii = ii(:)';
        for j = ii
            marker(j) = 1;
            qe = qe + 1;
            q(qe) = j;
        end
        qs = qs + 1;
    end
    compi = q(1:qe);
end
