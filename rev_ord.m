
function [iperm] = rev_ord(perm)
%% function [iperm] = rev_ord(perm)
     n = length(perm);
     for i=1:n
         iperm(perm(i)) = i;
     end
