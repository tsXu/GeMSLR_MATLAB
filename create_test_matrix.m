function [ A, n, gen, p, lev_ptr, subdm_ptr, nlev ] = create_test_matrix()
%[ A, n, gen, p, lev_ptr, subdm_ptr, nlev ] = create_test_matrix()
%   Create test matrix.
%   A:          return a test sparse matrix
%   n:          the size of the matrix
%   gen:        general matrix or not (for general matrix we don't have permutation in this test file yet)
%   p:          the permutation array
%   lev_ptr:    start of each level, plus the end of the last level + 1
%   subdm_ptr:  end of each subdomain on each level + 1

% choose to use general matrix or Laplacian
gen = false;
nlev = 2; % set number of levels
k = 4; % set number of subdomains on each level
minsep = max(k,16); % set minimal size of separator

if(gen)
    % load matrix data
    %load('./Matrices/young4c.mat');
    
    load('./Matrices/dwg961b.mat')
    %load('./Matrices/ted_AB.mat');
    
    % load permutation data
    % user need to load p, nlev, lev_ptr, and subdm_ptr
    %load('./Matrices/young4c_perm.mat');
    
    A = Problem.A;
    
    %A = mmread('./Matrices/bp__1000.mtx');
    n = size(A,1);
    [ p, nlev, lev_ptr, subdm_ptr ] = RKway_Gen( A, k, nlev, minsep );
    
else
    nx = 16;
    ny = 16;
    nz = 16;
    alphax = 0.0+0.0i;
    alphay = 0.0+0.0i;
    alphaz = 0.0+0.0i;
    shiftz = 0.0+0.0i;
    A = fd3d(nx,ny,nz,alphax,alphay,alphaz,shiftz);
    n = nx*ny*nz;
    %[p,~,~,lev_ptr,subdm_ptr] = NDRegGrid(nx, ny, nz, nlev); % grid
    [ p, nlev, lev_ptr, subdm_ptr ] = RKway_Gen( A, k, nlev, minsep ); % spectral
    
end

end

