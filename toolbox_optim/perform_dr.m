function [x,R] = perform_dr(x,ProxF,ProxG,options)

% perform_dr - Douglas Rachford algorithm
%
%   [x,R] = perform_dr(x,ProxF,ProxG,options);
%
%   Solves
%       min_x F(x)+G(x)
%   where F and G are proper convex functions, with an easy to compute
%   proximal mapping.
%
%   INPUTS:
%   ProxF(y,sigma) computes Prox_{sigma*F}(x)
%   ProxG(y,sigma) computes Prox_{sigma*G}(x)
%   options.niter is the number of iterations.
%   options.gamma the relaxation parameter for DR.
%   options.verb is for the diaplay of iterations.
%   options.report(x) is a function to fill in R.
%
%   OUTPUTS:
%   x is the final solution.
%   R(i) = options.report(x) at iteration i.
%
%   Copyright (c) 2010 Gabriel Peyre

rProxF = @(x,tau)2*ProxF(x,tau)-x;
rProxG = @(x,tau)2*ProxG(x,tau)-x;

mu = getoptions(options, 'mu', 1);
niter = getoptions(options, 'niter', 100);
verb = getoptions(options, 'verb', 1);
gamma = getoptions(options, 'gamma', 1);
report = getoptions(options, 'report', @(x)0);

clear R;
y = x;
for i=1:niter 
    % record energies
    if verb
        progressbar(i,niter);
    end
    R(i) = report(x);
        y = (1-mu/2)*y + mu/2*rProxF( rProxG(y,gamma),gamma );        
	x = ProxG(y,gamma);
end
