function [c,t] = load_curve(name, n, options)

% load_curve - load a sample curve
%
% function c = load_curve(name, n, options);
%
%   name can be 'circle' or 'square'
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
    otherwise
        error('Unkwnown curve');
end