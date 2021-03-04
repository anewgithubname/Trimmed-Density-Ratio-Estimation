function [l,g] = LLKLIEP(theta,kP,kQ,nu, lambda)
np = size(kP,2); nu_np = nu*np;
%sorting
t = theta'*kP; [~, id] = sort(t);

l = -sum(theta'*kP(:,id(1:nu_np)),2)./np + nu*log(mean(exp(theta'*kQ),2)) + lambda*theta'*theta;

N_q = sum(exp(theta'*kQ),2);
g_q = exp(theta'*kQ)./ N_q;
g = -sum(kP(:, id(1:nu_np)),2)./np + nu*kQ*g_q' + 2*lambda*theta;

end