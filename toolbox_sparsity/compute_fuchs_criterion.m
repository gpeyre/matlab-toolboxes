function [f,i,si] = compute_fuchs_criterion(D,x)

% compute_fuchs_criterion - compute the Fuchs criterion
%
%   D is an overcomplete dictionary of size (n,m), m>n
%   x is a m-dimensional vector
%
%   This criterion depends on sign(x) only
%
%   f = compute_fuchs_criterion(D,x);
%
%   f = max D(:,S1)' * pinv(D(:,S2)) * sign(x(S2))
%   
%   where   S2=find(abs(x)>0)
%           S1=find(abs(x)==0);
%
%   i the index (belongs to S1) where the maximum correlation is reached
%   si is the sign of the correlation.
%
%	f<1 => the vector can be recovered by L1 optimization
%
%	Inspired by 
%		J.J. Fuchs
%		Recovery of exact sparse representations in the presence of bounded noise.
%		IEEE-T-IT, vol.~51, 10, p.~3601--3608, oct.2005.
%
%   Important: the atoms of D must be unit-normed
%
%   Copyright (c) 2007 Charles Dossal and Gabriel Peyre

if size(x,2)>1 && size(x,1)>1    
    for k=1:size(x,2)
        [f(k),i(k),si(k)] = compute_fuchs_criterion(D,x(:,k));
    end
    return;
end

x = double( sign(x(:)) );
%% compute selected and un-selected columns
S1=find(abs(x)>0);
S2=find(x==0);
%% compute pseudo inverse
D1 = D(:,S1);
warning off;
D1 = (D1'*D1)^(-1) * D1';
warning on;
%% compute criterion
d0 = D(:,S2)' * D1' * x(S1);
[f,i] = max( abs( d0 ) );
si = sign(d0(i));
% convert to global indexes
i = S2(i); 
