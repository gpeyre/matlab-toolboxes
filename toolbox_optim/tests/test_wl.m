%%
% Add the toolbox.

addpath('../');
addpath('../toolbox/');

%% 
% Dimension of the signal, number of measurements and sparsity level.
n = 200;
m = n/2;
k = n/10;

%%
% Sensing matrix.
A = randn(m,n);

%%
% Sparse vector in [0,1].
done = 0;
while ~done
	I = randi(n,k,1);
	if length(unique(I)) == length(I), 
		done=1; 
	end
end
x0 = zeros(n,1);
x0(I) = 0.5*rand(size(I));
u0 = x0./(1-x0);

%%
% Measurements.
y = A*u0;


%%
% We want to solve 
% min_{A*u=y, u >= 0} ||u||_1 = min_{A*u=y, u >= 0} \sum_i u_i

%%
% This can be rewriten as the minimization of |F(x)+G(x)|
% where, F=\sum_i (u_i + i_{P_+}(u_i) and G=i_{A*x=y}.


%%
% The proximity operator of gamma F.
projPos = @(u)u.*(u>=0);
ProxF = @(u,gamma)projPos(u-gamma);

%%
% The proximity operator of the indicator of |A*x=y| is the orthogonal
% projection on A*u=y.
pA = A'*(A*A')^(-1);
ProxG = @(u,gamma)u + pA*(y-A*u);

%%
% Create a function to record the values of F and the constraint at each iteration.
F = @(u)sum(u);
Constr = @(u)1/2*norm(y-A*u)^2;
options.report = @(u)struct('F', F(u), 'Constr', Constr(u));

%%
% Run the algorithm. 
options.gamma = 50;
options.niter = 5000;
[u,R] = perform_dr(zeros(n,1), ProxF, ProxG, options);
u = projPos(u);
x = u./(1+u);
%uLP = linprog(ones(n,1),-eye(n,n),zeros(n,1),A,y);
%xLP = uLP./(1+uLP);

%%
% Display the solution. At convergence, it should be of sparsity |k|.
figure(1)
clf;
subplot(211)
stem(1:n,u,'.b');hold on
stem(I,u0(I),'or');hold off;axis tight;
legend('BP-DR','Original','Location','Best');legend boxoff
%stem([u uLP u0]);axis tight;
%legend('BP-DR','BP-IP','Original','Location','Best');legend boxoff
title('u=x/(1-x)');
subplot(212)
stem(1:n,x,'.b');hold on
stem(I,x0(I),'or');hold off;axis tight;
legend('BP-DR','Original','Location','Best');legend boxoff
%stem([x xLP x0]);axis tight;
%legend('BP-DR','BP-IP','Original','Location','Best');legend boxoff
title('x');

%%
% Retrieve the F and constraint function values.
f = s2v(R,'F');
constr = s2v(R,'Constr');

%%
% Display.
figure(2)
clf;
subplot(211);
plot(f(2:end));
axis tight; title('Objective');
subplot(212);
plot(constr(2:end));
axis tight; title('Constraint');
