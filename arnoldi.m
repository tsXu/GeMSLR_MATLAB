function [V, H, m2, tits] = arnoldi( Afunc, V, H, m, neig, neig2, orthtol, atol)
%  [V, H, m2, tits] = arnoldi( Afunc, V, H, m, neig, neig2, orthtol, atol)
%  Arnoldi process with re-orth and deflation
%  V:       is n x (m2) contains V (in and out), on entry V(:,1) is the
%  starting vector.
%  H:       is m2 x m2 contains H (in and out)
%  m2:      num of cols in V
%  tits:    total number of iterations
%  Afunc:   matrix of function that could compute matvec
%  m:       steps in the Arnoldi procedure
%  neig:    number of eigenvalues we want, stop when reached. Check every
%  10 loop
%  neig2:   number of eigenvalues we kept
%  orthtol: orthogonal tol
%  atol:    tol for eigenvalues
%%-----------------------------------------------

%% setup paras
n = size(V,1);
alpha   = 1/sqrt(2);
tits = 0;
check = -neig; % no need to check when its < neig
checked = false;

%main loop

V(:,1) = V(:,1)/norm(V(:,1),2);

for i = 1:m
    tits = tits + 1;
    i1 = i + 1;
    check = check + 1;
    
    % compute matvec
    V(:,i1) = feval(Afunc,V(:,i));

    % compute norm for re-orth
    normv = norm(V(:,i1));

    H(1:i,i) = 0;
    %MGS loop
    for j=1:i
        t = V(:,j)'*V(:,i1);
        H(j,i) = t;
        V(:,i1) = V(:,i1) - t*V(:,j);
    end
    
    % reorth
    t = norm(V(:,i1),2);
    while(t >= orthtol && t < alpha*normv)
        normv = t;
        %redo MGS loop
        for j=1:i
            t = V(:,j)'*V(:,i1);
            H(j,i) = H(j,i) + t;
            V(:,i1) = V(:,i1) - t*V(:,j);
        end
        t = norm(V(:,i1),2);
    end

    H(i1,i) = t;  
    if (t < orthtol)
        % if break, all eigenvalues are accurate
        fprintf("Lucky breakdown at step %d with %g\n",tits,t);
        break;
    end
    V(:,i1) = V(:,i1) / t;
    
    if (check >= 10)
        check = 0;
        %%-------------------  Compute decomposition
        m1 = tits;
        Hc = H(1:m1, 1:m1);
        Vc = V(:,1:m1);
        [X1,H1]  = schur(Hc);
        [X2,H2]  = eig(H1);
        X = X1*X2;
        vals = abs(diag(H2));
        
        %%-------------------  Compute Residual
        vres = H(m1+1,m1)*X(m1,:);
        
        
        %%-------------------- Select Convergenced Eigs
        cov_idx = find(abs(vres) < atol);
        
        if(length(cov_idx) >= neig)
            % we've got enough
            m2      = neig2;
            select = zeros(m1,1);
            [~, ord] = sort(vals(cov_idx),'descend');
            select(cov_idx(ord(1:m2)))=1;
            [XS,HS]  = ordschur(X1,H1,select);
            
            V       = zeros( n, m2);
            H       = zeros( m2, m2);
            
            H(1:m2,1:m2) = HS(1:m2,1:m2);
            V(:,1:m2) = Vc*XS(:,1:m2);
            
            checked = true;
            break;
        end
        
    end
    
end

if ~checked
    % haven't got enough
    m1 = tits;
    Hc = H(1:m1, 1:m1);
    Vc = V(:,1:m1);
    [X1,H1]  = schur(Hc);
    [X2,H2]  = eig(H1);
    X = X1*X2;
    vals = abs(diag(H2));
    
    %%-------------------  Compute Residual
    vres = H(m1+1,m1)*X(m1,:);
    
    %%-------------------- Select Convergenced Eigs
    cov_idx = find(abs(vres) < atol);
    ncov = length(cov_idx);
    
    if ncov == 0
        V = [];
        H = [];
        m2 = 0;
        return;
    end
    m2 = min(neig2, ncov);
    select = zeros(m1,1);
    [~, ord] = sort(vals(cov_idx),'descend');
    select(cov_idx(ord(1:m2)))=1;
    [XS,HS]  = ordschur(X1,H1,select);
    
    V       = zeros( n, m2);
    H       = zeros( m2, m2);

    H(1:m2,1:m2) = HS(1:m2,1:m2);
    V(:,1:m2) = Vc*XS(:,1:m2);
end

end
