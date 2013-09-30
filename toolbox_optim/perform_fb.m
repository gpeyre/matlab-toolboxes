function [x,R] = perform_fb(x, ProxF, GradG, L, options)

% perform_fb - Forward backward method
%
%    [x,R] = perform_fb(x, ProxF, GradG, L, options);
%
%   Solves min_x g(x) + f(x)
%   where g is a smooth convex proper function and f is a
%   convex proper function with an easy to compute proximal operator.
%
%   Use several first order-scheme depending on options.method:
%       options.method = 'fb' : classical Foward-backward
%       options.method = 'fista' : FISTA method of Beck and Teboule
%       options.method = 'nesterov' : Nesterov scheme.
%
%   INPUTS:
%   ProxF(y,sigma) computes Prox_{sigma*F}(x)
%   GradG(x) computes \nabla f(x)
%   L is the lipschitz constant of the gradient, if g is C^2:
%       L = max_x norm( Hg(x) ) 
%       where Hg(x) is the hessian of g at x. 
%       For instance, if g(x)=1/2*|A*x-y|^2 then tau = norm(A)^2.
%   options.niter is the number of iterations.
%   options.verb is for the diaplay of iterations.
%   options.report(x) is a function to fill in R.
%
%   OUTPUTS:
%   x is the final solution.
%   R(i) = options.report(x) at iteration i.
%
%   Copyright (c) 2010 Gabriel Peyre

options.null = 0;
method = getoptions(options, 'method', 'fb');
report = getoptions(options, 'report', @(x)0);
niter = getoptions(options, 'niter', 100);
verb = getoptions(options, 'verb', 1);
fbdamping = getoptions(options, 'fbdamping', 1.8);

clear R;
t = 1;  % fista & nesterov
tt = 2/L; gg = 0; A = 0; % nesterov
y = x;
x0 = x;
for i=1:niter 
  	R(i) = report(x);
    if verb
        progressbar(i,niter);
    end
    switch method
        case 'fb'
            x = ProxF( x-fbdamping/L*GradG(x), fbdamping/L );
        case 'fista'
            xnew = ProxF( y - 1/L*GradG(y), 1/L );
            tnew = (1+sqrt(1+4*t^2))/2;
            y = xnew + (t-1)/(tnew)*(xnew-x);
            x = xnew; t = tnew;
        case 'nesterov'
            a = (tt + sqrt(tt^2 + 4*tt*A))/2;
            v = ProxF( x0-gg, A );
            z = (A*x+a*v)/(A+a);
            x = ProxF( z - 1/L*GradG(z) , 1/L  );
            gg = gg +  a * GradG(x); % P'*(P*x-y);
            A = A + a;
        otherwise
            error('Unknown method');
            
    end      
end