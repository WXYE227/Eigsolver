function [Y,apply_times,niter_list,solve_runtime] = applyTrapFunc(R,X,c,shifts,gmres_tol)
m = length(c);
[n,col] = size(X);
Y = zeros(n,col);
RX = R(X);
[SRX,~,~,niter_list,apply_times,solve_runtime] = gmres_shift(R,RX,shifts,100,gmres_tol,100);
for it = 1:m
    Y = Y + c(it) * SRX(:,:,it);
end
end

