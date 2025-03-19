set(0,'defaulttextinterpreter','latex')
%%
load('d1.mat')
C_ = -C_;
e = eig(full(G_), full(C_));

figure(1)
scatter(real(e),imag(e))

figure(2)
subplot(1,2,1)
spy(G_)
subplot(1,2,2)
spy(C_)
%%
% number of eigenvalue inside
s = 20; 

center = -200+1000*1j;
r = 90;

opt = [];
opt.ncol = ceil(s*1.2);
opt.retol = 1e-8;
opt.itmax = 10;
opt.verbose = 0;
opt.gmrestol = 1e-9;
% as good as simple rule if grmestol = retol/1e2
opt.spurioustol = 1e-2;

opt.m = [8,8];
% one step converge
% opt.m = [8,32];

[Uout_c,lambda_out_c,solve_times_c,it_c,relerrs_c,solve_times_list_c,record_niter_list_c,solve_runtime_sum_c] = ...
trap_eigs_HSRR_sparse(G_,C_,center,r,opt);

opt.m = 64;

[Uout_s,lambda_out_s,solve_times_s,it_s,relerrs_s,solve_times_list_s,record_niter_list_s,solve_runtime_sum_s] = ...
trap_eigs_HSRR_sparse(G_,C_,center,r,opt);

% compare eig and composite solver
e_inside = e(abs(e-center)<r);
figure(3)
scatter(real(lambda_out_c),imag(lambda_out_c),'r*')
hold on 
scatter(real(e_inside),imag(e_inside),'bo')
axis equal
%%
figure(4)
semilogy(relerrs_c(1:it_c),'b-s','Markersize',10)
hold on
semilogy(relerrs_s(1:it_s),'r--*','Markersize',10)
hold off
legend('Composite','Simple')
ylabel('relative error')
xlabel('subspace iteration')
xticks([1,2])
title('$n_x=10$')
set(gca,'Fontsize',22)