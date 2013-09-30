function [x,R] = perform_dr_spingarn(x, ProxF, options)

% perform_dr_spingarn - Douglas Rachford algorithm for the sum of >=2 functions 
%
%   [x,R] = perform_dr_spingarn(x,ProxF,options);
%
%   Solves
%       min_x sum_i F_i(x)
%   where F_i are proper convex functions, with an easy to compute
%   proximal mapping.
%
%   INPUTS:
%   ProxF{i}(y,sigma) computes Prox_{sigma*F_i}(x)
%   options.niter is the number of iterations.
%   options.gamma the relaxation parameter for DR.
%   options.verb is for the diaplay of iterations.
%   options.report(x) is a function to fill in R.
%
%   OUTPUTS:
%   x is the final solution.
%   R(i) = options.report(x) at iteration i.
%
%   Copyright (c) 2011 Gabriel Peyre

m = length(ProxF);
s = size(x);
x = x(:);

%%
% initialize by copy
x0 = repmat(x, [1 m]);

report_old = getoptions(options, 'report', @(x)x);
options.report = @(x)report_old(reshape(x(:,1),s));

ProxF1 = @(x,gamma)ProxMult(x,gamma,ProxF);
ProxG  = @(x,gamma)repmat(mean(x,2), [1 m]);

[x,R] = perform_dr(x0,ProxF1,ProxG,options);
x = reshape(x(:,1),s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = ProxMult(x, gamma, ProxF)

for i=1:size(x,2)
    x(:,i) = ProxF{i}(x(:,i),gamma);
end

