function [k] = kernel_linear(x)
[d n] = size(x); 
tt = ones(d); tt = triu(tt,1); tt=tt(:); idx = tt~=0; 
k = zeros(sum(tt~=0),n);
for i = 1:n
    t = x(:,i) * x(:,i)';
    t = triu(t,1); t = t(:);
    k(:,i) = t(idx);
end
k = [k;x.^2];
end