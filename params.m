%%Create parameters for the test

%% ILU para
%ilu_setup.type = 'crout';
ilu_setup.type = 'ilutp';
ilu_setup.droptol = 1e-02;
ilu_setup.thresh = 0;

%% plot para
plotperm = true;
ploteig = false;
set(0, 'defaultaxesfontsize', 18);
FS = 'fontsize';    fs = 20;
FW = 'fontweight';  fw = 'bold';
LW = 'linewidth';   lw = 1.8;  lw2 = 2.8;
MS = 'markersize';  ms = 6;
IN = 'interpreter'; in = 'latex';
%image size
xx = 500;
yy = 320;

%% GMRES paras
kdim = 50;
tolits = 1e-06;
maxits = 600;

%% LOW-RANK paras
rank_k_top      = 20;   % target number of eigs computed on the top level
rank_k2_top     = 20;   % target number of eigs kept on the top level
rank_k_lower    = 20;  % target number of eigs computed on other levels
rank_k2_lower   = 20;  % target number of eigs kept on other levels
msteps_top      = 600;  % Arnoldi steps during each outer iteration on the top level
msteps_lower    = 600;  % Arnoldi steps during each outer iteration on other levels
atol_top        = 1e-12; % how accurate we compute those eigenvalues on the top level
atol_lower      = 1e-12; % how accurate we compute those eigenvalues on the other levels
% note: residual correction is applied to compute more accureat L and U
% solve