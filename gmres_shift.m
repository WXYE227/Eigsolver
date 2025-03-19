function [x,flag,relres,niter_list,apply,solve_runtime] = gmres_shift ...
    (Afun0,b,omega,restart,tol,maxiter)

% solves the systems (A-omega I)x = rhs
% ||(A-omega I)x - rhs||\(||rhs||) < e for all omega

% INPUTS
% A:     a function handle to apply a NxN matrix
% rhs:   a Nx1 vector
% omega: a vector containing the values of omega

% OUTPUTS
% x:       a NxNomega vector representing the solution, col(k) is for omega(k)
% iter:    The number of iterations needed, will be less than maxiter
% err:     ||(A-omega I)x - rhs||\(||n||), for each shift. err(k) is for omega(k)



if nargin < 6 || isempty(maxiter)
    maxiter = 100;
end

if nargin < 5 || isempty(tol)
    tol = 1e-6;
end

if nargin < 4 || isempty(restart)
    restart = 30;
end

if nargin < 3 || isempty(omega)
    omega = 0;
end

norms = @(X)sqrt(sum(abs(X).^2,1));


[N,Ncol] = size(b);
Nomega = length(omega);
Nrst = restart;
normb = norms(b);

if ~isa(Afun0,'function_handle')
    A = Afun0;
    Afun = @(x)A*reshape(x,N,[]);
else
    Afun = @(x)Afun0(reshape(x,N,[]));
end

solved = zeros(Ncol,1,Nomega);
solved_col = zeros(Ncol,1);
x = zeros(N,Ncol,Nomega);
relres = zeros(Ncol,1,Nomega);

Q = zeros(N,Nrst,Ncol);
H = zeros(Nrst+1,Nrst,Ncol);
vh10 = zeros(Nrst+1,Ncol);
r = zeros(N,1,Ncol);

x0 = zeros([size(b),Nomega]);
r0 = b;
h10 = norms(r0);
vh10(1,:) = h10;

apply = 0;

Q(:,1,:) = bsxfun(@rdivide,r0,h10);
indexes = 1:Ncol;

niter_list = zeros(1,Ncol);
solve_runtime = 0;
for j = 1:maxiter
    % generate krylov space
    k = j;
    index = indexes(solved_col == 0);
    apply = apply + length(index);
    tapply0 = tic;
    r(:,:,index) = Afun(Q(:,k,index));
    solve_runtime = solve_runtime + toc(tapply0);
    for ii = index
        Q_sub = reshape(Q(:,1:k,ii),[N,k]);
        r_sub = reshape(r(:,:,ii),[N,1]);
        H_sub = Q_sub' * r_sub;

        r_sub = r_sub - Q_sub * H_sub;
        temp = Q_sub' * r_sub;
        r_sub = r_sub - Q_sub * temp;
        H_sub = H_sub + temp;
        
        H(1:k,k,ii) = H_sub;
        r(:,:,ii) = r_sub;
        H(k+1,k,ii) = norm(r_sub);
    end
    % early stop
    for i = 1:length(index)
        ic = index(i);
        H_c = H(1:k+1,1:k,ic);
        for io = 1:Nomega
 
            if solved(ic,1,io) == 0
                H_omega = H_c - omega(io) * eye(k+1,k);
                y = H_omega \ vh10(1:(k+1),ic);
                z = vh10(1:(k+1),ic) - H_omega * y;
                relres(ic,1,io) = norm(z)/normb(ic);
                if relres(ic,1,io) < tol
                    solved(ic,1,io) = 1;
                    x(:,ic,io) = x0(:,ic,io) + reshape(Q(:,1:k,ic),N,[])*y;
                end
                
            end
            
        end
        if all(solved(ic,1,:))
            solved_col(ic) = 1;
            niter_list(ic) = j;
        end
    end
    
    if all(solved_col)
        flag = 0;
        return;
    end
    Q(:,k+1,index) = r(:,:,index)./H(k+1,k,index);
end
flag = 1;
end