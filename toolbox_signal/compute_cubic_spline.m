function y = compute_cubic_spline(x,n)

% compute_cubic_spline - evaluate the nth derivative of the cubic spline.
%
%   y = compute_cubic_spline(x,der);
%
%   der is the order of derivative.
%   set det to 0 (default)) for evaluation of the cubic spline function.
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<2
    n=0;
end

sx = sign(x);
x = abs(x) ;

I12 = (x>1)&(x<=2);
I01 = (x<=1);

if n==0
    y = I01.*( 2/3-x.^2.*(1-x/2) ) + I12.*( 1/6*(2-x).^3 );
elseif n==1
    y = I01.*( -2*x+3/2*x.^2 ) + I12.*( -1/2*(-2+x).^2 );
    y = y.*sx;
elseif n==2
    y = I01.*( -2+3*x ) + I12.*( 2-x );    
elseif n==3
    y = I01.*( 3 ) + I12.*( -1 );   
    y = y.*sx; 
else
    errror('Derivatives higher than n=3 not allowed.');
end