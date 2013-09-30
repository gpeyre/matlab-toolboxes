function fi = perform_rbf_interpolation(Xi,X,f)

% perform_rbf_interpolation - scattered data interpolation
%
%   fi = perform_rbf_interpolation(Xi,X,f);
%
%   It uses a linear spline radial basis function.
%       fi(x) = a*x + b + sum a_i * |x-x_i|
%   where x_i are the interpolation points
%   and a_i are computed to satisfy the interpolation property
%       fi(x_i)=f(x_i)
%
% X is the location of sampling location, of size (nbr_sampling_points,3)
% Xi is the location of interpolation location, of size (nbr_evaluation_points,3)
% f is the value if the function at the sampling location of size (nbr_sampling_points,1)
% fi is the interpolated values, of size (nbr_evaluation_points,1)
% 
% Code borrowed from Aranz 
% <http://aranz.com/research/modelling/theory/directmethods.html>
%
%   Copyright (c) 2005 Gabriel Peyré

if size(Xi,2)>size(Xi,1)
    Xi = Xi';
end
if size(X,2)>size(X,1)
    X = X';
end
f = f(:);

% find coefficient of the expansion
coeff = fitit(X,f);

% perform interpolation
for i=1:length(Xi)
    fi(i) = eval_direct( X, coeff, Xi(i,:) );
end

% Function which directly evaluates the RBF spline
%
% s(u)=linear poly +\sum_i=1^n coeff(i)*phi(u-cent(i,:) 
%
% at the point u where the linear polynomial is given  
% in terms of its monomial cofficients. 
%
% Syntax     s=eval_direct(cent,coeff,u)
%
% Inputs
%    cent   n by dim array of coordinates for the centres
%    coeff  n+dim +1 vector of coefficients with the linear
%           polynomial part last.
%    u      row vector    point at which to evaluate
%
% Output
%    s    Value of the RBF at position u
%
% Code is written for clarity rather than Matlab efficiency

function s = eval_direct(cent,coeff,u)

[n dim] = size(cent);
s=0;
for j=1:n
   s = s + coeff(j)*phi(  norm(u - cent(j,:) )   );
end
% Now the linear polynomial bit
s=s+coeff(n+1);               % The constant
for i=1:dim
   s=s+coeff(i+n+1)*u(i);     % The various components of u
end





% Function to find the radial basis function s consisting
% of a linear plus a sum of shifts of \phi( | x | )
% interpolating to data
% f_i at the point cent(i,:) for 1 \leq i \leq n.
%
% Syntax [coeff]=fitit(cent,f)
%
% Input
%
%    cent    n by dim array of centres
%    f       n by 1 vector of values at centres
%
% Output
%    coeff   (n +3) by 1 vector of coefficients with the 
%            coefficients of 1, x and y (2D) or 1, x, y 
%            and z (3D) last.
%
function [coeff]=fitit(cent,f)

A=direct(cent);
%
[n dim] =size(cent);
f=f(:);                % Make sure its a column
f=[f ;zeros(dim+1,1)]; % add zeros for the polynomial part
                       % at the end of the column
coeff=A\f;


% Form the (n+dim+1)*(n+dim+1) matrix corresponding to 
% RBF interpolation with basic function phi and linears.
%
% The centres are assumed given in the n by dim array cent.
% phi is assumed given as a function of r. It is coded in 
% the Matlab function phi. m
%
% Syntax  [A]=direct(cent)
%
% Input          
%         cent  n*dim array  coordinates of centers
%               of reals     
% Output  A     (n+dim+1)*   Symmetric matrix of the linear
%               (n+dim+1)    system obtained if we solve
%               array of     for the radial basis function
%                            interpolant directly.
%
% Write the matrix A in the form
%             B    P
%     A   =  
%             P^t  O
% where P is the polynomial bit.
%
function [A]=direct(cent)
[n dim]=size(cent);
A=zeros(n,n);
for i=1:n
    for j=1:i
        r=norm(cent(i,:)-cent(j,:));
        temp=phi(r);
        A(i,j)=temp;
        A(j,i)=temp;
    end
end
%
% Now the polynomial part
%
P=[ones(n,1) cent];
A = [ A      P
      P' zeros(dim+1,dim+1)];


  

%
% Function phi of r
%
% Syntax [u] = phi(r)
%
% Remember if using something like the thinplate spline
 
% in 2D you will need to test for r nonzero before 
% taking the log.
function u=phi(r)
u =r;