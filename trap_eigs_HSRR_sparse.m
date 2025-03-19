function [Uout,lambda_out,solve_times,it,relerrs,solve_times_list,record_niter_list,solve_runtime_sum] = trap_eigs_HSRR_sparse(A,B,center,r,varargin)
% using sparese LU
solve_times = 0;
[ncol,m,reltol,itmax,verbose,gmres_tol,spurioustol,lusetting] = trap_options(varargin{:});
n_useful = ncol;
emax = abs(center) + r;
n = size(A,1);
if m(1) ~= 0
    [poles,weights] = trap_coef(center,r,m(1));
    prefact_L = cell(m(1),1);
    prefact_U = cell(m(1),1);
    prefact_p = cell(m(1),1);
    prefact_q = cell(m(1),1);
    t0=tic;
    for i = 1:m(1)
        [prefact_L{i},prefact_U{i},prefact_p{i},prefact_q{i}] = lu(poles(i)*B-A,lusetting,'vector');
    end
    toc(t0)
    R = @(X)applyTrapMat_sparse(prefact_L,prefact_U,prefact_p,prefact_q,B,weights,X,A,poles,verbose);
    
    compose = (length(m) == 2)
    if compose
        [u_poles,~] = trap_coef(0,1,m(2));
        c = u_poles./(1 + u_poles) / m(2);
        shifts = 1./(1 + u_poles);
        RR = @(X)applyTrapFunc(R,X,c,shifts,gmres_tol);
    else
        RR = R;
    end
else
    disp('m1 must not equal to 0')
end

% examine the eigenvalue of precodition_A and gmres_shift
if compose && verbose
    figure(10)
    disp('---------start eig pre A----------')
    pre_A = R(eye(n,n));
    [~,D] = eig(pre_A);
    lambda = diag(D);
    [poles,~] = trap_coef(0,1,m(2));
    gmres_poles = 1 ./ (1 + poles);
    
    scatter(real(lambda),imag(lambda),'b')
    hold on
    scatter(real(gmres_poles),imag(gmres_poles),'r')
    legend('precondition\_A','gmres\_shifts')
    hold off
    disp('---------end eig pre A----------')
end
%
rng(0)
Y = randn(n,ncol) + randn(n,ncol) * 1j;
[Y,~] = qr(Y,0);

% use for record the ture eigenvalue of 
idx = ones(ncol,1);
filtered_num_ = ncol; 
relerrs = ones(1,itmax);
solve_times_list = zeros(1,itmax);
record_niter_list = zeros(itmax,ncol);
solve_runtime_sum = 0;
for it = 1:itmax
    % disp(it)
    idx = ones(n_useful,1);
    % apply rational filter to V
    if compose
        [U,apply,niter_list,solve_runtime] = RR(Y);
        solve_runtime_sum = solve_runtime_sum + solve_runtime;
        % disp(max(niter_list))
        record_niter_list(it,1:n_useful) = niter_list;
    else
        ttmp = tic;
        [U,apply] = RR(Y);
        solve_runtime_sum = solve_runtime_sum + toc(ttmp);
    end

    % the total iter times should be niter * n_useful
    if compose
        solve_times_list(it) = apply * m(1);
    else
        solve_times_list(it) = max(m(1),1) * n_useful;
    end
    solve_times = solve_times + solve_times_list(it);
    
    % construct and solve the subprolblem
    [V,R] = qr(U,0);
    absR = abs(diag(R));
    max_absR = max(absR);
    
    % examine if it is sigular
    for i = 1:n_useful
        if absR(i)/max_absR < 1e-15
            n_useful = n_useful - 1;
            idx(i) = 0;
        end
    end
    
    V = V(:,idx==1);
    [W,~] = qr(A*V-center*(B*V),0);
    
    A_m = W' * A * V;
    B_m = W' * B * V; 
    [H_A,H_B,~,P_R] = qz(A_m,B_m);
    [V_R,D] = eig(H_A,H_B);
    lambda = diag(D);
    evs = V * (P_R * V_R);

%     temp = lsqminnorm(B*V,A*V);
%     norm(A*V-B*V*temp,'fro')/norm(A*V,'fro')
    % index the true eigenvalues, the eigenvalue outside of have huge error
    AX = A * evs; BX = B * evs; LBX = BX .* (lambda.');
    residual = vecnorm(AX-LBX) ./ (vecnorm(BX)*emax);
%     idx = (residual < spurioustol);
    idx = (residual < spurioustol) .* (abs(lambda - center) < r).';
    disp(['not spuriousnum',num2str(sum(residual < spurioustol))])
    % estimate for the convergence of true eigenvalues
    filtered_num = nnz(idx);
    
    residual_inside = residual(idx==1);
    lambda_inside = lambda(idx==1);
    [~, sort_index] = sort(abs(lambda_inside-center));
    % residual_sort = residual_inside(sort_index)
    if filtered_num ~= 0
        relerrs(it) = max(residual(idx==1));
        % estimate if the filtered_num is fixed
        if (filtered_num == filtered_num_ || it == 1) && relerrs(it)<reltol
            break;
        end
    end
    filtered_num_ = filtered_num;
    Y = evs;
end
Uout = evs(:,idx==1);
lambda_out = lambda(idx==1);
% solve_runtime_sum
end

