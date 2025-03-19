function [Y, niter] = applyTrapMat_sparse(prefact_L,prefact_U,prefact_p,prefact_q,B,weights,X,A,poles,verbose)
m = length(prefact_L);
BX = B * X;
TX = zeros(size(X));
Y = zeros(size(X));
e = zeros(1,m);
for i = 1:m
    S = prefact_U{i}\ (prefact_L{i}\BX(prefact_p{i},:));
    Y(prefact_q{i},:) = Y(prefact_q{i},:) + weights(i) * S;
    
    if verbose
        TX(prefact_q{i},:) = S;
        ABX = (poles(i)*B - A) * TX;
        e(i) = norm(ABX - BX,'fro')/norm(BX(prefact_p{i},:),'fro');
    end
end
if verbose
    disp(e)
end
niter = 0;
end