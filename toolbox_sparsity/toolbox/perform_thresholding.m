function y = perform_thresholding(x, t, type)

% perform_thresholding - perform hard or soft thresholding
%
%   y = perform_thresholding(x, t, type);
%
%   type is either 'hard' or 'soft' or 'semisoft'
%   t is the threshold
%
%   works also for complex data, and for cell arrays.
%
%   if type is 'strict' then it keeps the t largest entry in each
%   column of x.
%
%   Copyright (c) 2006 Gabriel Peyre

if nargin<3
    type = 'hard';
end

if iscell(x)
    % for cell arrays
    for i=1:size(x,1)
        for j=1:size(x,2)
            y{i,j} = perform_thresholding(x{i,j},t, type);
        end
    end
    return;
end

switch lower(type)
    case {'hard', ''}
        y = perform_hard_thresholding(x,t);
    case 'soft'
        y = perform_soft_thresholding(x,t);
    case 'semisoft'
        y = perform_semisoft_thresholding(x,t);
    case 'strict'
        y = perform_strict_thresholding(x,t);
    otherwise
        error('Unkwnown thresholding type.');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = perform_strict_thresholding(X,s)
%% keep only the s largest coefficients in each column of X

v = sort(abs(X)); v = v(end:-1:1,:);
v = v(round(s),:);
v = repmat(v, [size(X,1) 1]);
X = X .* (abs(X)>=v);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = perform_hard_thresholding(x,t)

t = t(1);
y = x .* (abs(x) > t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = perform_soft_thresholding(x,t)

if not(isreal(x))
    % complex threshold
    d = abs(x);
    d(d<eps) = 1;
    x = x./d .* perform_soft_thresholding(d,t);
end

t = t(1);
s = abs(x) - t;
s = (s + abs(s))/2;
y = sign(x).*s;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = perform_semisoft_thresholding(x,t)

if length(t)==1
    t = [t 2*t];
end
t = sort(t);

y = x;
y(abs(x)<t(1)) = 0;
I = find(abs(x)>=t(1) & abs(x)<t(2));
y( I ) = sign(x(I)) .* t(2)/(t(2)-t(1)) .* (abs(x(I))-t(1));