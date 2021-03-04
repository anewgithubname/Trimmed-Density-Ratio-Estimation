clear; rng(15);
d = 20; n = 500;

rng(23); adj = rand(d)<.1; adj = triu(adj, 1) + triu(adj, 1)'; thetaP = double(adj)*.45 + eye(d)*2;
rng(23); adj = rand(d)<.1-8e-2; adj = triu(adj, 1) + triu(adj, 1)'; thetaQ = double(adj)*.45 + eye(d)*2;
xp = mvnrnd(zeros(1,d),inv(thetaP),n)'; xp = [xp, [10*ones(d/2,1);-10*ones(d/2,1)]];
xq = mvnrnd(zeros(1,d),inv(thetaQ),n)';
figure; imagesc(thetaP-thetaQ)

kp = kernel_linear(xp); kq = kernel_linear(xq);

nu = (size(xp,2)-10)/size(xp,2)
% nu = 1;

theta = sparse(zeros(size(kq,1),1));
lambda = .7*log(d)/sqrt(n);
tic
step = 1; slength = inf; iter = 0; fold = inf;
while(slength > 1e-5)
    [f, gt] = LLKLIEP(theta,kp,kq, nu, 0);
    g = zeros(size(gt));
    
    id = abs(theta)>0;
    g(id) = gt(id) + lambda*sign(theta(id));
    id = theta==0 & gt > lambda;
    g(id) = gt(id) - lambda;
    id = theta==0 & gt < -lambda;
    g(id) = gt(id) + lambda;
    theta = theta - step*g./(iter+1);
    slength = step*norm(g)./(iter+1);
    fdiff = abs(f - fold);
    
    theta(abs(theta)<1e-4) = 0;
    theta = sparse(theta);
    %display some stuffs
    if iter > 50000
        disp('max iteration reached.')
        break;
    else
        iter = iter+1;
        fdiff = abs(f - fold);
        fold = f;
        if ~mod(iter,100)
            disp(sprintf('%d, %.5f, %.5f, %.5f, nz: %d',...
                iter, slength,fdiff,full(fold),full(sum(theta(1:end-d)~=0))))
        end
    end
end
toc

r_hat = @(theta, kt, kq) exp(theta'*kt - log(mean(exp(theta'*kq))));
rp = r_hat(theta, kp, kq);


%%
tt = ones(d); tt = triu(tt,1); tt=tt(:); idx = tt~=0; edges = [];
Delta = zeros(1,d*d);
Delta(idx) = abs(theta(1:end-d)); Delta = reshape(Delta,d,d);
Delta = Delta + Delta';
hfig = figure; hfig.PaperPosition = [0 0 8.5 8.5*.7776]; hfig.Position = [488 333 642 500]
imagesc(Delta); colorbar; h = title(sprintf('Estimated KL %.3f',-f)); hold on;
h.FontSize = 18;

GT = (thetaP - thetaQ);
for i = 1:d
    for j = 1:d
        if GT(i,j)~= 0
            h = rectangle('Position',[i-.5,j-.5,1,1]); h.EdgeColor = 'r'; h.LineWidth = 4
        end
    end
end
h = gca; h.FontSize = 18; axis normal; h.GridAlpha = .5;
print('change','-dpng')
