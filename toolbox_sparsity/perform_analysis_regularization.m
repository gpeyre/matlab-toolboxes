function [g,g_list,E] = perform_analysis_regularization(f, G, options)

% perform_analysis_regularization - perform a sparse regularization
%
%   [g,g_list,E] = perform_analysis_regularization(f, A, options);
%
% Method solves, given f of length n, for
%       min_g  E(g) = 1/2*|f-g|^2 + lambda * |A*g|_1
% where A is an (m,n) operator and |x|_1=\sum_k |x(k)|
%
%   A can be a (m,n) matrix or an operator
%       y = A(x,dir,options);
%   where dir=1 to computer y=A*x and y=A'*x when dir=-1.
%
%   You have to set lambda in options.lambda.
%   The number of iteration is options.niter.
%   The regularization parameter is options.eta (should be small enough).
%   You can use an increasing lambda by setting options.lambda_min and
%   options.lambda_max. In that case, g_list contains the list of solutions
%   for increasing values of lambda.
%
%   It uses the projection algorithm described in 
%       A. Chambolle
%       "An algorithm for Total variation Minimization and applications"
%       Journal of Mathematical Imaging and Vision, 20(1), p. 89-97, 2004
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
if isfield(options, 'niter')
    niter = options.niter;
else
    niter = 100;
end
if isfield(options, 'eta')
    eta = options.eta;
else
    eta = 0.1;
end
if isfield(options, 'g')
    h = options.h;
else
    h = f*0;
end

if isfield(options, 'lambda')
    lambda = options.lambda;
else
    lambda = .3;
end
if isfield(options, 'lambda_min')
    lambda_min = options.lambda_min;
else
    lambda_min = lambda;
end
if isfield(options, 'lambda_max')
    lambda_max = options.lambda_max;
else
    lambda_max = lambda;
end
if isfield(options, 'verb')
    verb = options.verb;
else
    verb = 1;
end

lambda_list = linspace(lambda_min,lambda_max,niter);

f = f(:);
n = length(f);

if nargout>2
    g_list = zeros(n,niter);
end

E = [];
p = appG(G,h,options);
for i=1:niter
    if verb
        progressbar(i,niter);
    end
    lambda = lambda_list(i);
    h = lambda*appGT(G,p,options);
    pp = -1/lambda * appG(G,h-f,options);
    % update
    p = ( p+eta*pp )./( 1+eta*abs(pp) );
    if nargout>2
        g = f-h;
        pg = appG(G,g,options);
        E(end+1) = .5*norm(h, 'fro')^2 + lambda_max * sum( abs(pg(:)) );
    end
    if nargout>1
        g_list(:,i) = f - h;
    end
end
g = f - h;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=appG(G,x, options)
if isnumeric(G)
    y = G*x;
else
	y = feval( G,x, 1, options );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=appGT(G,x, options)
if isnumeric(G)
    y = G'*x;
else
	y = feval( G,x, -1, options );
end
