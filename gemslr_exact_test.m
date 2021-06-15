%%This is a test file for GeMSLR
% Note: this is a version for complex linear system
% see params.m to change settings
% see create_test_matrix.m to modify test matrices
close all
clear

% read test parameters
params;

use_scinv = false;

% reset random seed
rng(0);

% create test matrix

[ A, n, gen, p, lev_ptr, subdm_ptr, nlev ] = create_test_matrix();

% create initial guess and rhs
rhs = A * ones(n,1);
x0 = sin(1:n)';

% apply permutation
A0 = A;
rhs0 = rhs;
A = A(p,p);
rhs = rhs(p);

% plot permutation matrix or not
if(plotperm)
    hFig=figure(1); spy(A0);
    xlabel('');
    title(sprintf('Original matrix'),IN,in);
    set(hFig, 'Position', [500 0 xx yy])
    hFig=figure(2); spy(A); hold on;
    for i=1:nlev-1
        subdmi_ptr = subdm_ptr{i};
        subdmi_ptr = subdmi_ptr(:)';
        for j=subdmi_ptr
            vert_hori_lines(n,j,'--','color', [0.0 0.5 0.0],'LineWidth',1.5);
        end
    end
    for i=2:nlev
        vert_hori_lines(n,lev_ptr(i),'r','LineWidth',2);
    end
    xlabel('');
    set(hFig, 'Position', [1000 0 xx yy])
    Problem.name = 'Single-level vertex separator';
    title(sprintf('%s', Problem.name),IN,in);
end

fprintf(1, '----preprocessing done....\n');

%% GeMSLR, ILUT phase

% compute exact factorization
ilu_setup.droptol = 0e-02;

% preset fill factor
PREnnzLU = 0;
PREnnzLR = 0;

tic;
PRE.nlev = nlev;
for levi = 1:nlev-1
    % on each level
    % | B_i   F_i |
    % | E_i   C_i |
    
    % Extract B_i, F_i, and E_i
    n_s = lev_ptr(levi);
    n_e = lev_ptr(levi+1)-1;
    PRE.Levs{levi}.A  = A(n_s:n,n_s:n);
    PRE.Levs{levi}.B  = A(n_s:n_e,n_s:n_e);
    PRE.Levs{levi}.F  = A(n_s:n_e, n_e+1:n);
    PRE.Levs{levi}.E  = A(n_e+1:n, n_s:n_e);
    PRE.Levs{levi}.C  = A(n_e+1:n, n_e+1:n);
    
    % ILU of B
    [PRE.Levs{levi}.LB, PRE.Levs{levi}.UB] = ilu( PRE.Levs{levi}.B, ilu_setup);
    PREnnzLU = PREnnzLU + nnz(PRE.Levs{levi}.LB) + nnz(PRE.Levs{levi}.UB);
end

% ILUT of C on the last level
n_s = lev_ptr(nlev);
PRE.Levs{nlev}.A  = A(n_s:n, n_s:n);
[ PRE.Levs{nlev}.LB, PRE.Levs{nlev}.UB] = ilu(PRE.Levs{nlev}.A, ilu_setup);
PREnnzLU = PREnnzLU + nnz(PRE.Levs{nlev}.LB) + nnz(PRE.Levs{nlev}.UB) - n;
toc;
fprintf('ILU time\n');

%% Build of low-rank correction
tic;
for levi = (nlev-1):-1:1
    % # of nodes on current level
    nB = size(PRE.Levs{levi}.LB, 1);
    % # of nodes on remaining levels
    nC = n - lev_ptr(levi+1) + 1;
    % set rank for this level
    if levi == 1
        rank_k  = rank_k_top; % compute
        rank_k2 = min(rank_k,rank_k2_top); % keep
        msteps  = msteps_top;
        atol    = atol_top;
    else
        rank_k  = rank_k_lower; % compute
        rank_k2 = min(rank_k,rank_k2_lower); % keep
        msteps  = msteps_lower;
        atol    = atol_lower;
    end
    
    if (msteps > 0 && rank_k2 > 0)
        
        % can't compute more eigs than we have
        rank_ki = min(rank_k, nC);
        mstepsi = min(msteps, nC);
        
        % apply standard Arnoldi
        Arnodi_V = zeros(nC,rank_ki+1);
        Arnodi_V(:,1) = randn(nC,1);
        %Arnodi_V(:,1) = ones(nC,1);
        Arnodi_H = zeros( rank_ki+1,rank_ki+1);
        
        if use_scinv
            I = eye(nC);
            S = PRE.Levs{levi}.C - PRE.Levs{levi}.E * (PRE.Levs{levi}.B \ PRE.Levs{levi}.F);
            SCinv = I - (PRE.Levs{levi}.C - PRE.Levs{levi}.E * ( PRE.Levs{levi}.B \ PRE.Levs{levi}.F ) ) / PRE.Levs{levi}.C;
            Afunc = @(x) axpy( SCinv, x);
        else
            SB = PRE.Levs{levi}.B - PRE.Levs{levi}.F * (PRE.Levs{levi}.C \ PRE.Levs{levi}.E);
            CESF = PRE.Levs{levi}.C \ (PRE.Levs{levi}.E * ( SB \ PRE.Levs{levi}.F ) );
            Afunc = @(x) axpy( CESF, x);
        end
        
        [Arnodi_V, Arnodi_H, m2, tits] = arnoldi( Afunc, Arnodi_V, Arnodi_H, msteps, rank_k, rank_k2, 1e-16, atol);
        fprintf('SCinv Find %d eigs out of %d on lev %d with %d iterations\n',m2 ,nC ,levi, tits);
        
        if m2 == 0
            % can't find any eigenvalue, return with empty low-rank correction
            PRE.Levs{levi}.G = zeros(0,0);
            PRE.Levs{levi}.Z = zeros(nC,0);
            PRE.Levs{levi}.Y = zeros(nC,0);
        else
            % the rank we keep in this approximation
            local_rank = m2;
            
            % map to the low-rank correction
            theta = 0.0;
            if use_scinv
                PRE.Levs{levi}.G = inv(eye(local_rank) - Arnodi_H) - 1/(1-theta)*eye(local_rank);
            else
                PRE.Levs{levi}.G = Arnodi_H;
            end
            % here we keep Z and Y for future tests,
            % we are not counting them twice when compute the fill factor
            % when you change them remember to change the fill factor
            PRE.Levs{levi}.Z = Arnodi_V;
            PRE.Levs{levi}.Y = Arnodi_V;
            
            fprintf('Rank number is %d on lev %d\n',local_rank ,levi);
        end
        
    else
        % no need to build correction
        PRE.Levs{levi}.G = zeros(0,0);
        PRE.Levs{levi}.Z = zeros(nC,0);
        PRE.Levs{levi}.Y = zeros(nC,0);
    end
    
    % Z is same as Y, only count once
    PREnnzLR = PREnnzLR + nnz(PRE.Levs{levi}.Z) + nnz(PRE.Levs{levi}.G);
end
toc;
fprintf('Build low-rank time\n');

% preconditioner solve function
PrecFunc = @(x) solve_levi(PRE,x,1);
PrecFunc2 = @(x) solve_levi_no_low_rank(PRE,x,1);

% get fill-factor
PREnnz = PREnnzLU + PREnnzLR;
ffact1 = PREnnz/nnz(A);
ffact2 = PREnnzLU/nnz(A);
%% Complex GMRES
fprintf('Starting iterations with low-rank...\n');
tic;
n = size(A,1);
[sol1,res1,its1] = fgmrez(A,PrecFunc,rhs,x0,tolits,maxits,kdim);
toc;

%% Complex GMRES without low-rank
fprintf('Starting iterations without low-rank...\n');
tic;
%n = size(A,1);
[sol2,res2,its2] = fgmrez(A,PrecFunc2,rhs,x0,tolits,maxits,kdim);
toc;

%% output final result
% plot convergence result
figure;
semilogy([0:its1],res1,'linestyle','--','marker','*','LineWidth',2,'color','r')
text(its1,res1(its1)/2,num2str(ffact1),'fontsize',18);
title([' SLR of matrix -- N = ', num2str(n)]);
hold on;
semilogy([0:its2],res2,'linestyle','--','marker','*','LineWidth',2,'color','b')
text(its2,res2(its2)/2,num2str(ffact2),'fontsize',18);

legend('GeMSLR-fgmres', 'NoLowRank-fgmres');
% print message
fprintf('GEMSLR...\n');
fprintf('Its: %d, ||b-Ax||/||b|| = %e\n', its1, norm(rhs-A*sol1)/norm(rhs));
fprintf('Fill-in ratio is %f (lu=%f, lrk=%f)\n', PREnnz/(nnz(A)), PREnnzLU/nnz(A), PREnnzLR/nnz(A));
fprintf('GEMSLR without low-rank update...\n');
fprintf('Its: %d, ||b-Ax||/||b|| = %e\n', its2, norm(rhs-A*sol2)/norm(rhs));

%% update solution
x(p) = sol1;
