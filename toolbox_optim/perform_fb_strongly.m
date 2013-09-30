function [x,R] = perform_fb_strongly(x, K, KS, GradFS, ProxGS, L, options)

% perform_admm - preconditionned ADMM method
%
%    [x,R] = perform_fb_strongly(x, K, KS, GradFS, ProxGS, L, options);
%
%   Minimization of 
%       min_x F(x)+G(K*x)                   (*)
%   where F is strongly convex, and F is a proper convex function.
%
%   Use the equivalence of (*) with
%       min_u F^*(-K^* u) + G^*(u)          (**)
%   using x = grad(F^*)(-K^* u)
%
%   Where the convex dual function is
%       F^*(y) = sup_x <x,y>-F(x)
%
%   Uses perform_fb to solve (**).
%
%   INPUTS:
%   GradFS(x) is grad(F^*)(x)
%   ProxGS(u,tau) is prox_{tau*G^*}(u)
%   K(x) is a linear operator
%   KS(u) is K^* the dual of K.
%   L is the lipshitz constant of K*GradFS*KS
%   options.niter is the number of iterations.
%   options.verb is for the diaplay of iterations.
%   options.report(x) is a function to fill in R.
%
%   OUTPUTS:
%   x=grad(F^*)(-K^* u) is the final solution.
%   R(i) = options.report(x) at iteration i.
%
%   Copyright (c) 2010 Gabriel Peyre

if isnumeric(K)
    K = @(x)K*x;
end  
if isnumeric(KS)
    KS = @(x)KS*x;
end   

report_old = getoptions(options, 'report', @(x)0);
options.report = @(u)report_old(GradFS(-KS(u)));

NewGrad = @(u)-K(GradFS(-KS(u)));
u0 = K(x);
[u, R] = perform_fb(u0, ProxGS, NewGrad, L, options);
x = GradFS(-KS(u));
