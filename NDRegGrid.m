function [p,t,tlev,lev_ptr,subdm_ptr] = NDRegGrid(nx, ny, nz, lev)
n = nx * ny * nz;
A = fd3d(nx, ny, nz, 0, 0, 0, 0);
if (nz < 2)
    [p,t,tlev] = NDRegGrid2D([nx,ny],lev,1,[1:n]',[],[]);
else
    [p,t,tlev] = NDRegGrid3D(A, [nx,ny,nz],lev,1,[1:n]',[],[]);
end
% p has the ''interleaving'' ND order
t = [1;t];
nt = length(t);
for i=2:nt
    t(i) = t(i) + t(i-1);
end
t = [t(1:nt-1), t(2:nt)-1];
% now each row of t is a row range of each tree node
% sort the tree nodes in t by levels (high to low)
[tlev,idx] = sort(tlev, 'descend');
t = t(idx,:);
% the desired ND order
p2 = [];
lev_ptr = 0;
tlev = [tlev; 0];
ii = 1;
subdmi_ptr = [];
for i=1:nt-1
    p2 = [p2; p(t(i,1):t(i,2))];
    subdmi_ptr = [subdmi_ptr, length(p2)+1];
    if (tlev(i) ~= tlev(i+1))
        lev_ptr = [lev_ptr; length(p2)];
        subdm_ptr{ii} = subdmi_ptr;
        subdmi_ptr = [];
        ii = ii + 1;
    end
end
lev_ptr = lev_ptr + 1;
p = p2;
end

function [p,t,tlev] = NDRegGrid2D(size, lev, ilev, p, t, tlev)
if (ilev < lev)
    [p1,p2,p3,size1,size2] = RegGrid2DPart(size, p);
    [p1,t,tlev] = NDRegGrid2D(size1, lev, ilev+1, p1, t, tlev);
    [p2,t,tlev] = NDRegGrid2D(size2, lev, ilev+1, p2, t, tlev);
    p = [p1; p2; p3];
    node_size = length(p3);
else
    node_size = size(1)*size(2);
end
t = [t; node_size];
tlev = [tlev; ilev];
end

function [p1,p2,p3,size1,size2] = RegGrid2DPart(size, p0)
nx = size(1);
ny = size(2);
n = nx*ny;
if (nx <= ny)
    % partition along X
    assert(ny >= 3);
    cut = round(ny/2);
    ny1 = cut - 1;
    ny2 = ny - cut;
    %
    n1 = ny1 * nx;
    n2 = ny2 * nx;
    n3 = nx;
    %
    par1 = [1:n1];
    par2 = [n-n2+1:n];
    par3 = [n1+1:n1+n3];
    %
    size1 = [nx, ny1];
    size2 = [nx, ny2];
else
    % partition along Y
    assert(nx >= 3);
    cut = round(nx/2);
    nx1 = cut - 1;
    nx2 = nx - cut;
    %
    n1 = nx1 * ny;
    n2 = nx2 * ny;
    n3 = ny;
    %
    id = [0:n-1]';
    idx = floor(id/ny);
    idy = mod(id, ny);
    perm = idy * nx + idx + 1;
    %
    par1 = perm(1:n1);
    par2 = perm(n-n2+1:n);
    par3 = perm(n1+1:n1+n3);
    %
    size1 = [ny, nx1];
    size2 = [ny, nx2];
end
p1 = p0(par1); 
p2 = p0(par2); 
p3 = p0(par3);
end



%% 3-D
function [p,t,tlev] = NDRegGrid3D(A, size, lev, ilev, p, t, tlev)
if (ilev < lev)
    [p1,p2,p3,size1,size2] = RegGrid3DPart(size, p);
    A3 = A(p3,p3);  q = amd(A3); p3 = p3(q);
    [p1,t,tlev] = NDRegGrid3D(A, size1, lev, ilev+1, p1, t, tlev);
    [p2,t,tlev] = NDRegGrid3D(A, size2, lev, ilev+1, p2, t, tlev);
    p = [p1; p2; p3];
    node_size = length(p3);
else
    AA = A(p,p);  q = amd(AA);  p = p(q);
    node_size = size(1)*size(2)*size(3);
end
t = [t; node_size];
tlev = [tlev; ilev];
end

function [p1,p2,p3,size1,size2] = RegGrid3DPart(size, p0)
nx = size(1);
ny = size(2);
nz = size(3);
nxy = nx*ny;
nxz = nx*nz;
nyz = ny*nz;
n = nx*ny*nz;

if (nxy <= nxz && nxy <= nyz)
    % partition along XY
    assert(nz >= 3);
    cut = round(nz/2);
    nz1 = cut - 1;
    nz2 = nz - cut;
    %
    n1 = nxy * nz1;
    n2 = nxy * nz2;
    n3 = nxy;
    %
    par1 = [1:n1];
    par2 = [n-n2+1:n];
    par3 = [n1+1:n1+n3];
    %
    size1 = [nx, ny, nz1];
    size2 = [nx, ny, nz2];
elseif (nxz <= nxy && nxz <= nyz)
    % partition along XZ
    assert(ny >= 3);
    cut = round(ny/2);
    ny1 = cut - 1;
    ny2 = ny - cut;
    %
    n1 = nxz * ny1;
    n2 = nxz * ny2;
    n3 = nxz;
    % id: x z y order
    id = [0:n-1]';
    idy = floor(id/nxz);
    idxz = mod(id, nxz);
    idz = floor(idxz/nx);
    idx = mod(idxz, nx);
    perm = idz*nxy + idy*nx + idx + 1;
    %
    par1 = perm(1:n1);
    par2 = perm(n-n2+1:n);
    par3 = perm(n1+1:n1+n3);
    %
    size1 = [nx, nz, ny1];
    size2 = [nx, nz, ny2];
elseif (nyz <= nxy && nyz <= nxz)
    % partition along YZ
    assert(nx >= 3);
    cut = round(nx/2);
    nx1 = cut - 1;
    nx2 = nx - cut;
    %
    n1 = nyz * nx1;
    n2 = nyz * nx2;
    n3 = nyz;
    % id: y z x order
    id = [0:n-1]';
    idx = floor(id/nyz);
    idyz = mod(id, nyz);
    idz = floor(idyz/ny);
    idy = mod(idyz, ny);
    perm = idz*nxy + idy*nx + idx + 1;
    %
    par1 = perm(1:n1);
    par2 = perm(n-n2+1:n);
    par3 = perm(n1+1:n1+n3);
    %
    size1 = [ny, nz, nx1];
    size2 = [ny, nz, nx2];
end

p1 = p0(par1); 
p2 = p0(par2); 
p3 = p0(par3);
end