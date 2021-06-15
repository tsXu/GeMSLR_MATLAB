function [ map ] = spectralkway( A, k )
%[ map ] = spectralkway( A, k )
%   partition the adjancency graph of |A|+|A'| into k parts.
%   in this simple test we assume that A in connected
%   k must be some power of 2
%   A: matrix
%   k: parts
%   map: array of domains each node belons to. Starts from 1.

if floor(k)~=k || k < 1
    fprintf("Input k must be positive integer.\n");
    map = [];
    return;
end

% k = 2^k1
k1 = log(k)/log(2);
if abs(round(k1)-k1) > 1e-06
    fprintf("Input k must be power of 2.\n");
    map = [];
    return;
end
k1 = round(k1);

if k1 == 0
    n = size(A,1);
    map = ones(1,n);
    return;
end

%% get graph Laplacian
[ ~, W, ~ ] = unweighted_laplacian( A );

%% start partition
map = lapkway( W, k1, 1);

end

function [map] = lapkway( W, tlvl, clvl)
    
    if size(W,1) == 0
        map = [];
    end

    % create laplacian
    n = size(W,1);
    D = sum(W,2);
    D = spdiags(D,0,n,n);
    L = D - W;
    
    map = zeros(n,1);
    
    % get all the connect components
    [comps] = connected_components(L);
    ncomp = length(comps);
    
    for compi = 1:ncomp
        
        ni = comps{compi};
        nni = length(ni);
        Li = L(ni, ni);
        Wi = W(ni, ni);
        
        % get the second smallest eigenvector
        opts.p = min(nni, 20);
        
        if (nni > 1)
            try
                if(nni > 32)
                    [v,d] = eigs(Li, 2, 'sm', opts);
                    [~,i] = max(abs(diag(d)));
                    v = v(:,i);
                else
                    [v,d] = eig(full(Li));
                    [~,i] = max(abs(diag(d)));
                    v = v(:,i);
                end
            catch
                warning('Problem using function.  Assigning a value of 0.');
                v = randn(nni, 1);
            end

            % try to equallly split
            vmed = median(v);
            imed = find(v == vmed);
            nmed = length(imed);
            randi = randn(nmed, 1);
            v(imed(randi>0)) = vmed + 1;
            v(imed(randi<=0)) = vmed - 1;
            
            i1 = find(v>vmed);
            i2 = find(~(v>vmed));
        else
            i1 = 1:nni;
            i2 = [];
        end
        
        W1 = Wi(i1,i1);
        W2 = Wi(i2,i2);
        
        if clvl < tlvl
            if ~isempty(i1)
                map(ni(i1)) = lapkway( W1, tlvl, clvl+1);
            end
            
            if ~isempty(i2)
                map(ni(i2)) = lapkway( W2, tlvl, clvl+1) + (tlvl - clvl) * 2;
            end
        else
            % last level
            map(ni(i1))=1;
            map(ni(i2))=2;
        end
    end
end

