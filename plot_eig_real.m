function [] = plot_eig_real(vals)
%plot_eig_real(vals)
%   vals: eigenvalues

n = length(vals);

figure;
hold on;

for i = 1:n
    vali = real(vals(i));
    plot([vali,vali],[-1,1],'b');
end

set(gcf, 'Position',  [0, 0, 800, 100])

end

