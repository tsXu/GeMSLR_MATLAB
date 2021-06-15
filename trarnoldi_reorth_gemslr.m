function [W,R,m2,tits] = trarnoldi_reorth_gemslr (PRE, levi, v, m, maxits, neig, tol, resits)
%  [W,R,m2,res] = arnoldi_reorth_gemslr (PRE, levi, v, m, maxits, neig, tol)
%  arnoldi process with re-orth, thick-restart, and deflation
%  try to find eigs clost to val
%  W:       is n x (m2)  contains convergenced eigs
%  R:       is m2 x m2
%  m2:      is the real number of eigs we found
%  tits:    total number of iterations
%  PRE:     preconditioning data structure
%  levi:    level to solve
%  v:       initial guess
%  maxits:  max number of iterations
%  neig:    number of eigenvalues we want
%  tol:     convergence tol
%  resits:  number of residual correction when solve with B and C
%%-----------------------------------------------

%% setup paras
n       = size(v,1);
orthtol = min(1e-16,tol/100);
%alpha1  = 1e-02;
mm      = min(m+neig,n);
V       = zeros(n, mm+1);
H       = zeros(mm+1,mm);
V(:,1)  = v / norm(v,2);

tr_factor = 0.15;
num_cov = 0;
its     = 0;
tsteps  = 0;

trlen   = 0;

tits = 0;
while its < maxits
    its = its + 1;
    %fprintf('Its number %d\n',its);
    k = trlen + 1;
    
    %% arnoldi start from k
    msteps   = min(m+num_cov,n);
    %[its, msteps]
    %[maxits,msteps]
    [V,H,m2,tits1] = arnoldi_standard_defl(PRE, levi, V, H, k, msteps, [], [], 0, orthtol, resits);
    tits = tits + tits1;
    tsteps = tsteps + m2 - k;
    %spy(H);
    %pause();
    %r = A*V(:,1:m2) - V(:,1:m2)*H(1:m2,1:m2) - V(:,m2+1) * H(m2+1,1:m2);
    %fprintf('Arnoldi residual norm %e\n',norm(r));
    %sum(r,1)
    %norm(r)
    
    %% restart
    %restart detect
    %%-------------------- Arnoldi Objects 
    Vm = V(:,1:m2);
    Hm = H(1:m2,1:m2);
    %%-------------------- Get Ritz values and vectors
    if num_cov > 0
        %Hmu = Hm(1:num_cov,1:num_cov);
        %Hml = Hm(num_cov+1:m2,num_cov+1:m2);
        %H12 = Hm(1:num_cov,num_cov+1:m2);
        %O12 = zeros(num_cov,m2-num_cov);
        %O21 = zeros(m2-num_cov,num_cov);
        %[Xml,Hml]  = schur(Hml);
        %Iml = eye(num_cov);
        %X1 = [Iml, O12; O21, Xml];
        %H1 = [Hmu, H12*Xml; O21, Hml];
        [X1,H1]  = schur(Hm);
    else
        [X1,H1]  = schur(Hm);
    end
    [X2,H2]  = eig(H1);
    X = X1*X2;
    vals = diag(H2);
    vecs = Vm*X;
    %%-------------------  Compute Residual
    vres = H(m2+1,m2)*X(m2,:);
    
    %%-------------------- selece convergenced eigs
    cov_idx = find(abs(vres) < tol);
    %if(length(cov_idx)<num_cov)
    %    without "lock", number of convergenced eigenvalues might reduce
    %    during the iteration
    %    warning('reduce');
    %end
    num_cov = length(cov_idx);
    
    ico_idx = find(abs(vres) >= tol);
    num_ico = length(ico_idx);
    
    if (num_cov >= neig || its >= maxits)
        %fprintf('Arnoldi break at its %d of %d, find %d eigenvalues\n', its, maxits, ncov);
        break;
    end
    
    %% restart
    if (num_ico == 0)
        V      = zeros(n,mm+1);
        H      = zeros(mm+1,mm);
        
        V(:,1) = rand(n,1);
        
        V(:,1) = V(:,1)/norm(V(:,1));
        trlen = 0;
    elseif (num_ico == 1)
        %%-------------------- restart with the only option
        V      = zeros(n,mm+1);
        H      = zeros(mm+1,mm);
        v      = vecs(:,ico_idx);
        
        V(:,1) = Vm*v;
        
        V(:,1) = V(:,1) / norm(V(:,1));
        trlen = 0;
    else    
        %%-------------------- thick restart
        %trlen = min(5,num_ico);
        nicov_pick = min(ceil(num_ico*tr_factor),num_ico);
        npick = max(1,nicov_pick);
        %trlen = 2+num_cov;
        trlen = npick+num_cov;
        select = zeros(m2,1);
        
        %%-------------------- thick restart with the left and right most
        %[~, ord] = sort(vals(ico_idx));
        %select(ico_idx(ord(1))) = 1;
        %select(ico_idx(ord(num_ico))) = 1;
        
        %or other options
        [~, ord] = sort(vals(ico_idx),'descend');
        select(ico_idx(ord(1:npick))) = 1;
        select(cov_idx)=2;
        [XS,HS]  = ordschur(X1,H1,select);
        vlast  = V(:,m2+1);
        beta   = H(m2+1,m2);
        V      = zeros(n,mm+1);
        H      = zeros(mm+1,mm);
        H(1:trlen,1:trlen) = HS(1:trlen,1:trlen);
        
        V(:,1:trlen) = Vm*XS(:,1:trlen);
        
        nvlast = norm(vlast);
        V(:,trlen+1) = vlast/nvlast;
        H(trlen+1,num_cov+1:trlen) = beta * nvlast * XS(m2,num_cov+1:trlen);
        %r = A*V(:,1:trlen) - V(:,1:trlen)*H(1:trlen,1:trlen) - vlast * H(trlen+1,1:trlen);
        %fprintf('Restart ncov %d, residual norm %e\n',num_cov,norm(r));
        %'d'
    end
    
    %% restart over
    
end

select = zeros(m2,1);
select(cov_idx)=1;
[XS,~]  = ordschur(X1,H1,select);
m2=num_cov;
W = Vm*XS(:,1:num_cov);
R = zeros(m2,m2);
for i = 1:num_cov
    %R(1:i,i) = W(:,1:i)'*(EBinvFCinv(PRE,W(:,i),levi,resits));
    R(1:i,i) = H1(1:i,i);
end

end
