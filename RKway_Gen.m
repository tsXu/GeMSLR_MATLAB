function [ p, nlev, lev_ptr, subdm_ptr ] = RKway_Gen( A, k, nlev, minsep )
%[ p, nlev, lev_ptr, subdm_ptr ] = RKway_Gen( A, nlev, minsep )
%   Apply recursive K-way partition (spectral-based). We assume that the
%   adjancency graph of A has only 1 connected component.
%   return:
%   p:          the permutation array
%   nlev:       number of actual levels
%   lev_ptr:    start of each level, plus the end of the last level + 1
%   subdm_ptr:  end of each subdomain on each level + 1
%   inputs:
%   A:          target matrix
%   k:          the k of Kway partition
%   nlev:       target number of levels
%   minsep:     min separator (stop when matrix size is smaller)

n = size(A,1);
clev = 1;

B = A;

% initialize unit permutation
p = 1:n;
lev_ptr = 1;
n_shift = 0;

% we can't partition matrix to more than its num of rows
minsep = max(minsep, k);

while( clev < nlev && n > minsep)
    clev = clev + 1;
    lev_ptr = [lev_ptr, lev_ptr(end)];
    subdm_ptr{clev-1}=zeros(k,1);
    % local permutation array
    lp = zeros(n,1);
    idx_s = 1;
    idx_e = n;
    
    % spectral partition
    [ map ] = spectralkway( B, k);
    
    % call METIS partition
    %map = metiskway(B,k);
    
    % loop through each components
    for i = 1:k
        % get nodes
        nodes = find(map==i);
        nnodes = map~=i;
        lennodes = length(nodes);
        for j = 1:lennodes
            row = nodes(j);
            if(nnz(B(row,nnodes)))
                % not empty, exterior node
                lp(idx_e) = row + n_shift;
                idx_e = idx_e - 1;
            else
                % empty, interior node
                lev_ptr(clev) = lev_ptr(clev) + 1;
                lp(idx_s) = row + n_shift;
                idx_s = idx_s + 1;
            end
        end
        subdm_ptr{clev-1}(i) = lev_ptr(clev);
    end
    
    % ready for the next level
    p(n_shift+1:end) = p(lp);
    n_shift = lev_ptr(clev) - 1;
    B=A(p(n_shift+1:end),p(n_shift+1:end));
    n = size(B,1);
end

nlev = clev;
n = size(A,1);
subdm_ptr{clev} = n+1;
lev_ptr = [lev_ptr, n+1];

for i = 1:clev
    % apply RCM on each subdomain
    for j = 1:length(subdm_ptr{i})
        % get diagonal matrix
        if j == 1
            idx_s = lev_ptr(i);
            idx_e = subdm_ptr{i}(1)-1;
            C = A(p(idx_s:idx_e),p(idx_s:idx_e));
            rcm = symrcm(C);
            p(idx_s:idx_e) = p(rcm + idx_s - 1);
        else
            idx_s = subdm_ptr{i}(j-1);
            idx_e = subdm_ptr{i}(j)-1;
            C = A(p(idx_s:idx_e),p(idx_s:idx_e));
            rcm = symrcm(C);
            p(idx_s:idx_e) = p(rcm + idx_s - 1);
        end
        
    end
end

end

