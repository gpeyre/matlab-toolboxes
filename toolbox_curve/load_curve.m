function [c,t] = load_curve(name, n, options)

% load_curve - load a sample curve
%
% c = load_curve(name, n, options);
%
%   name can be 'circle', 'square', 'spiral', 'cloudy', 'rand', 'levy'.
%
%   some parameters such as regularity (alpha) can be modified through
%   options.alpha=xxx.
%
%   Copyright (c) 2004 Gabriel Peyré

options.null = 0;

if isfield(options, 'sampling')
    sampling = options.sampling;
else
    sampling = 'uniform';
end

switch lower(sampling)
    case 'uniform'
        t = linspace( 0,1-1/n, n );
        t = t(:);
    case 'random'
        t = sort( rand(n,1) );
    case 'parabolic'
        t = linspace( -1,1, n+1 ); t = t(1:end-1)';
        t = (1+t.^3)/2;
    otherwise
        error('Unkwnown sampling');
end

c = zeros(n,2);
switch lower(name)
    case 'circle'
        t = t*2*pi;
        c = [cos(t), sin(t)];
    case 'square'
        I = find(t<0.25);
        c(I,:) = [ t(I)*4, t(I)*0 ];
        I = find(t>= 0.25 & t<0.5);
        c(I,:) = [ t(I)*0+1 (t(I)-0.25)*4 ];
        I = find(t>= 0.5 & t<0.75);
        c(I,:) = [ 1-(t(I)-0.5)*4, t(I)*0+1 ];
        I = find(t>= 0.75);
        c(I,:) = [ t(I)*0, 1-(t(I)-0.75)*4 ];
    case 'spiral'
        tr = 3; % nbr tours
        r = log(tr*t+1);
        r = tr*t+1;
        c = [cos(2*pi*tr*t).*r, sin(2*pi*tr*t).*r];
    case 'cloudy'
        c = exp( 2i*pi*t ) + 0.05*exp( 2i*pi*10*t )  + 0.0*exp( 2i*pi*20*t );
        c = [ real(c), imag(c)];
    case 'rand'
        if isfield(options, 'alpha')
            alpha = options.alpha;
        else
            alpha = 1.5;
        end
        c = gen_rand_curve(n, alpha, options);
    case 'levy'
        type = 'isotropic';
        if isfield(options, 'alpha')
            alpha = options.alpha;
        else
            alpha = 1.2;
        end
        if isfield(options, 'sigma')
            sigma = options.sigma;
        else
            sigma = 0.1;
        end
        beta = 0; delta = 0; % no drift
        c = gen_levy_flight(n,alpha,sigma,beta,delta,type);
    otherwise
        error('Unkwnown curve');
end





function c = gen_rand_curve(n, alpha, options)

% gen_rand_curve - generate a curve with C^alpha smoothness
%
%   c = gen_rand_curve(n, alpha, options);
%
%   Copyright (c) 2005 Gabriel Peyré

n = round(n/2)*2;

options.null = 0;
if isfield(options, 'rmax')
    rmax = options.rmax;
else
    rmax = 0.4;
end

% amplitud deacy
h = (-n/2+1):(n/2);
h = (abs(h)+1).^(-alpha-0.5);
h = rescale( h, 0,rmax );
% add random phase
y = h .* exp( rand(1,n)*2i*pi );

y = fftshift(y);
y(1) = 0;   % zero mean
% around a circle
y(2) = 1;
y(n) = 0;

y = ifft(n*y);

% close the curve
y(end+1) = y(1);

c = [real(y); imag(y)];