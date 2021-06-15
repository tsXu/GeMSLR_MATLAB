function [] = vert_hori_lines(n,jj,varargin)
% [] = vert_hori_lines(n,jj,varargin)
% plot two lines
j = jj-0.5;
plot([1,n],[j,j],varargin{:});
plot([j,j],[1,n],varargin{:});
end
