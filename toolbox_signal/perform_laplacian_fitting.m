function [sigma,alpha,oor,oor1,oor2] = perform_laplacian_fitting( h, name, options )

% perform_laplacian_fitting - compute the parameter for a laplcian distribution
%
%   [a,b] = perform_laplacian_fitting( h, name, otions );
%
%   name can be :
%       'genlaplacian' : exp(-|x/sigma|^alpha)
%       'laplacian' : exp(-|x/sigma|)
%       'genlaplacian0' : only positive entries
%       'laplacian0' : only positive entries
%
%   Copyright (c) Gabriel Peyré 2006

options.null = 0;
if isfield(options, 'alpha_precision')
    a = options.alpha_precision;
else
    a = 64;
end
if isfield(options, 'sigma_precision')
    s = options.sigma_precision;
else
    s = 64;
end

m = length(h);

if isfield(options, 'alpha_min')
    alpha_min = options.alpha_min;
else
    alpha_min = 0.2;
end
if isfield(options, 'alpha_max')
    alpha_max = options.alpha_max;
else
    alpha_max = 2;
end

if isfield(options, 'sigma_min')
    sigma_min = options.sigma_min;
else
    sigma_min = 0.05;
end
if isfield(options, 'sigma_max')
    sigma_max = options.sigma_max;
else
    sigma_max = 0.2;
end

if strcmp(name, 'laplacian') || strcmp(name, 'laplacian0')
    alpha_min = 1; alpha_max = 1; a = 1;
end

sigma = linspace(sigma_min,sigma_max,s);
alpha = linspace(alpha_min,alpha_max,a);

if isfield(options, 'sigma')
    sigma = options.sigma(:)'; s = length(sigma);
end
if isfield(options, 'alpha')
    alpha = options.alpha(:)'; a = length(alpha);
end

% Initial estimate of alpha by moment matching
if 0
    if name(end)=='0'
        u = linspace(0,1,length(h))';
    else
        u = linspace(-1,1,length(h))';
    end
    m1 = sum( h(:) .* u );
    m2 = sum( h(:) .* u.^2 );
    alpha1 = estbeta(m1, m2);    
    
    [tmp,ialpha] = min( abs(alpha-alpha1) );
    alpha = alpha(ialpha);
    a = 1;
end

H = compute_laplacian_distribution( name, m, sigma, alpha );

% compute best fit for Culback distance
%  -log2(H) .* h
I = find(H<eps); H(I)=eps;
C = sum( -log2(H) .* repmat(h,[1,s,a]), 1 );
[tmp,i] = min(C(:));
[isigma,ialpha] = ind2sub( [s,a], i );

sigma = sigma(isigma);
alpha = alpha(ialpha);

oor1 = (isigma==1 || isigma==s) && m>2;
oor2 = (ialpha==1 || ialpha==a) && a>1 && m>2;
oor = oor1 || oor2;
    
return;

if (isigma==1 || isigma==s) && m>2
    warning('Maybe out of range fit.');
end
if (ialpha==1 || ialpha==a) && a>1 && m>2
    warning('Maybe out of range fit.');
end






function beta = estbeta(m1, m2)
% ESTBETA Estimate the beta parameter of a generalized Gaussian
%       beta = estbeta(m1, m2)
%
% Input:
%	m1: derivation from the mean
%	m2: variance of input data
%
% Output:
%	beta: paramaters of generalized Gaussian p.d.f.
%	      f(x) = K * exp(-(abs(x) / alpha)^beta)
%
% Used by GGMLE.  See also GGMME
%
% Note:	This is a faster version that computes beta using interpolation
%	where SBPDF use FZERO to solve inverse function
%
% Author: Minh N. Do, Dec. 1999

% beta = F^(-1)(m1^2 / m2);

x = 0.05 * [1:100]';
f = [2.4663e-05 0.004612 0.026349 0.062937 0.10606 0.15013 0.19234 0.23155 ...
     0.26742 0.3 0.32951 0.35624 0.38047 0.40247 0.42251 0.44079 0.45753 ...
     0.47287 0.48699 0.5 0.51202 0.52316 0.53349 0.5431 0.55206 0.56042 ...
     0.56823 0.57556 0.58243 0.58888 0.59495 0.60068 0.60607 0.61118 ...
     0.616 0.62057 0.6249 0.62901 0.63291 0.63662 0.64015 0.64352 ...
     0.64672 0.64978 0.65271 0.6555 0.65817 0.66073 0.66317 0.66552 ...
     0.66777 0.66993 0.672 0.674 0.67591 0.67775 0.67953 0.68123 0.68288 ...
     0.68446 0.68599 0.68747 0.68889 0.69026 0.69159 0.69287 0.69412 ...
     0.69532 0.69648 0.6976 0.69869 0.69974 0.70076 0.70175 0.70271 ...
     0.70365 0.70455 0.70543 0.70628 0.70711 0.70791 0.70869 0.70945 ...
     0.71019 0.71091 0.71161 0.71229 0.71295 0.71359 0.71422 0.71483 ...
     0.71543 0.71601 0.71657 0.71712 0.71766 0.71819 0.7187 0.7192 0.71968];

val = m1^2 / m2;
if val < f(1)
    beta = 0.05;
elseif val > f(end)
    beta = 5;
else
    beta = interp1q(f, x, val);
end