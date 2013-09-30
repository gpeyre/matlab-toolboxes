function y = perform_cpx_dualtree_transform(x,Jmin,options)

% perform_cpx_dualtree_transform - perform a complex valued dual-tree wavelet transform
%
%   y = perform_cpx_dualtree_transform(x,Jmin,options);
%
%   Copyright (c) 2004 Gabriel Peyré

options.null = 0;

if iscell(x)
    dir = -1;
    n = size(x{1},1)*2;
else
    dir = 1;
    n = size(x,1);
end

Jmax = log2(n)-1;
J = Jmax-Jmin+1;

[Faf, Fsf] = FSfarras;
[af, sf] = dualfilt1;

if dir==1
    y1 = cplxdual2D(x, J, Faf, af);
    % flatten everything
    y = {};
    for j=1:J
        for i=1:2
        for d1=1:2
            for d2=1:3
                y{end+1} = y1{j}{i}{d1}{d2};
            end
        end
        end
    end
    for d1=1:2
        for i=1:2
            y{end+1} = y1{end}{i}{d1};
        end
    end
else
    % flatten everything
    x1 = {};
    k = 0;
    for j=1:J
        for i=1:2
        for d1=1:2
            for d2=1:3
                k = k+1;
                x1{j}{i}{d1}{d2} = x{k};
            end
        end
        end
    end
    for d1=1:2
        for i=1:2
        k = k+1;
        x1{J+1}{i}{d1} = x{k};
        end
    end
    y = icplxdual2D(x1, J, Fsf, sf);
end


function w = cplxdual2D(x, J, Faf, af)

% Dual-Tree Complex 2D Discrete Wavelet Transform
%
% USAGE:
%   w = cplxdual2D(x, J, Faf, af)
% INPUT:
%   x - 2-D array
%   J - number of stages
%   Faf{i}: first stage filters for tree i
%   af{i}:  filters for remaining stages on tree i
% OUTPUT:
%   w{j}{i}{d1}{d2} - wavelet coefficients
%       j = 1..J (scale)
%       i = 1 (real part); i = 2 (imag part)
%       d1 = 1,2; d2 = 1,2,3 (orientations)
%   w{J+1}{m}{n} - lowpass coefficients
%       d1 = 1,2; d2 = 1,2 
% EXAMPLE:
%   x = rand(256);
%   J = 5;
%   [Faf, Fsf] = FSfarras;
%   [af, sf] = dualfilt1;
%   w = cplxdual2D(x, J, Faf, af);
%   y = icplxdual2D(w, J, Fsf, sf);
%   err = x - y;
%   max(max(abs(err)))
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

% normalization
x = x/2;

for m = 1:2
    for n = 1:2
        [lo w{1}{m}{n}] = afb2D(x, Faf{m}, Faf{n});
        for j = 2:J
            [lo w{j}{m}{n}] = afb2D(lo, af{m}, af{n});
        end
        w{J+1}{m}{n} = lo;
    end
end

for j = 1:J
    for m = 1:3
        [w{j}{1}{1}{m} w{j}{2}{2}{m}] = pm(w{j}{1}{1}{m},w{j}{2}{2}{m});
        [w{j}{1}{2}{m} w{j}{2}{1}{m}] = pm(w{j}{1}{2}{m},w{j}{2}{1}{m});
    end
end


function y = icplxdual2D(w, J, Fsf, sf)

% Inverse Dual-Tree Complex 2D Discrete Wavelet Transform
% 
% USAGE:
%   y = icplxdual2D(w, J, Fsf, sf)
% INPUT:
%   w - wavelet coefficients
%   J - number of stages
%   Fsf - synthesis filters for final stage
%   sf - synthesis filters for preceeding stages
% OUTPUT:
%   y - output array
% See cplxdual2D
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

for j = 1:J
    for m = 1:3
        [w{j}{1}{1}{m} w{j}{2}{2}{m}] = pm(w{j}{1}{1}{m},w{j}{2}{2}{m});
        [w{j}{1}{2}{m} w{j}{2}{1}{m}] = pm(w{j}{1}{2}{m},w{j}{2}{1}{m});
    end
end

y = zeros(size(w{1}{1}{1}{1})*2);
for m = 1:2
    for n = 1:2
        lo = w{J+1}{m}{n};
        for j = J:-1:2
            lo = sfb2D(lo, w{j}{m}{n}, sf{m}, sf{n});
        end
        lo = sfb2D(lo, w{1}{m}{n}, Fsf{m}, Fsf{n});
        y = y + lo;
    end
end

% normalization
y = y/2;




function y = sfb2D(lo, hi, sf1, sf2)

% 2D Synthesis Filter Bank
%
% USAGE:
%   y = sfb2D(lo, hi, sf1, sf2);
% INPUT:
%   lo, hi - lowpass, highpass subbands
%   sf1 - synthesis filters for the columns
%   sf2 - synthesis filters for the rows
% OUTPUT:
%   y - output array
% See afb2D
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/


if nargin < 4
    sf2 = sf1;
end

% filter along rows
lo = sfb2D_A(lo,    hi{1}, sf2, 2);
hi = sfb2D_A(hi{2}, hi{3}, sf2, 2);

% filter along columns
y = sfb2D_A(lo, hi, sf1, 1);


function [af, sf] = FSfarras

% Farras filters organized for the dual-tree
% complex DWT.
%
% USAGE:
%    [af, sf] = FSfarras
% OUTPUT:
%    af{i}, i = 1,2 - analysis filters for tree i
%    sf{i}, i = 1,2 - synthesis filters for tree i
% See farras, dualtree, dualfilt1.
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

af{1} = [
                  0                  0
  -0.08838834764832  -0.01122679215254
   0.08838834764832   0.01122679215254
   0.69587998903400   0.08838834764832
   0.69587998903400   0.08838834764832
   0.08838834764832  -0.69587998903400
  -0.08838834764832   0.69587998903400
   0.01122679215254  -0.08838834764832
   0.01122679215254  -0.08838834764832
                  0                  0
 ];
   
sf{1} = af{1}(end:-1:1, :);

af{2} = [
   0.01122679215254                  0
   0.01122679215254                  0
  -0.08838834764832  -0.08838834764832
   0.08838834764832  -0.08838834764832
   0.69587998903400   0.69587998903400
   0.69587998903400  -0.69587998903400
   0.08838834764832   0.08838834764832
  -0.08838834764832   0.08838834764832
                  0   0.01122679215254
                  0  -0.01122679215254
];

sf{2} = af{2}(end:-1:1, :);


function [af, sf] = dualfilt1

% Kingsbury Q-filters for the dual-tree complex DWT
%
% USAGE:
%    [af, sf] = dualfilt1
% OUTPUT:
%    af{i}, i = 1,2 - analysis filters for tree i
%    sf{i}, i = 1,2 - synthesis filters for tree i
%    note: af{2} is the reverse of af{1}
% REFERENCE:
%    N. G. Kingsbury,  "A dual-tree complex wavelet
%    transform with improved orthogonality and symmetry
%    properties", Proceedings of the IEEE Int. Conf. on
%    Image Proc. (ICIP), 2000
% See dualtree
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

% These cofficients are rounded to 8 decimal places.

af{1} = [
   0.03516384000000                  0
                  0                  0
  -0.08832942000000  -0.11430184000000
   0.23389032000000                  0
   0.76027237000000   0.58751830000000
   0.58751830000000  -0.76027237000000
                  0   0.23389032000000
  -0.11430184000000   0.08832942000000
                  0                  0
                  0  -0.03516384000000
 ];
 
af{2} = [
                  0  -0.03516384000000
                  0                  0
  -0.11430184000000   0.08832942000000
                  0   0.23389032000000
   0.58751830000000  -0.76027237000000
   0.76027237000000   0.58751830000000
   0.23389032000000                  0
  -0.08832942000000  -0.11430184000000
                  0                  0
   0.03516384000000                  0
];
 
sf{1} = af{1}(end:-1:1, :);
 
sf{2} = af{2}(end:-1:1, :);



function [lo, hi] = afb2D(x, af1, af2)

% 2D Analysis Filter Bank
%
% USAGE:
%   [lo, hi] = afb2D(x, af1, af2);
% INPUT:
%   x - N by M matrix
%       1) M, N are both even
%       2) M >= 2*length(af1)
%       3) N >= 2*length(af2)
%   af1 - analysis filters for columns
%   af2 - analysis filters for rows
% OUTPUT:
%    lo - lowpass subband
%    hi{1} - 'lohi' subband
%    hi{2} - 'hilo' subband
%    hi{3} - 'hihi' subband
% EXAMPLE:
%   x = rand(32,64);
%   [af, sf] = farras;
%   [lo, hi] = afb2D(x, af, af);
%   y = sfb2D(lo, hi, sf, sf);
%   err = x - y;
%   max(max(abs(err)))
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

if nargin < 3
   af2 = af1;
end

% filter along columns
[L, H] = afb2D_A(x, af1, 1);

% filter along rows
[lo,    hi{1}] = afb2D_A(L, af2, 2);
[hi{2}, hi{3}] = afb2D_A(H, af2, 2);



function [lo, hi] = afb2D_A(x, af, d)

% 2D Analysis Filter Bank
% (along one dimension only)
%
% [lo, hi] = afb2D_A(x, af, d);
% INPUT:
%    x - NxM matrix, where min(N,M) > 2*length(filter)
%           (N, M are even)
%    af - analysis filter for the columns
%    af(:, 1) - lowpass filter
%    af(:, 2) - highpass filter
%    d - dimension of filtering (d = 1 or 2)
% OUTPUT:
%     lo, hi - lowpass, highpass subbands
%
% % Example
% x = rand(32,64);
% [af, sf] = farras;
% [lo, hi] = afb2D_A(x, af, 1);
% y = sfb2D_A(lo, hi, sf, 1);
% err = x - y;
% max(max(abs(err)))

lpf = af(:, 1);     % lowpass filter
hpf = af(:, 2);     % highpass filter

if d == 2
   x = x';
end

N = size(x,1);
L = size(af,1)/2;
x = cshift2D(x,-L);

lo = upfirdn(x, lpf, 1, 2);
lo(1:L, :) = lo(1:L, :) + lo([1:L]+N/2, :);
lo = lo(1:N/2, :);

hi = upfirdn(x, hpf, 1, 2);
hi(1:L, :) = hi(1:L, :) + hi([1:L]+N/2, :);
hi = hi(1:N/2, :);

if d == 2
   lo = lo';
   hi = hi';
end




function y = cshift2D(x, m)

% 2D Circular Shift
% 
% USAGE:
%    y = cshift2D(x, m)
% INPUT:
%    x - M by N array
%    m - amount of shift
% OUTPUT:
%    y - matrix x will be shifed by m samples down
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

[N, M] = size(x);
n = 0:N-1;
n = mod(n-m, N);
y = x(n+1,:);



function y = sfb2D_A(lo, hi, sf, d)

% 2D Synthesis Filter Bank
% (along single dimension only)
%
% y = sfb2D_A(lo, hi, sf, d);
% sf - synthesis filters
% d  - dimension of filtering
% see afb2D_A


lpf = sf(:, 1);     % lowpass filter
hpf = sf(:, 2);     % highpass filter

if d == 2
   lo = lo';
   hi = hi';
end

N = 2*size(lo,1);
L = length(sf);
y = upfirdn(lo, lpf, 2, 1) + upfirdn(hi, hpf, 2, 1);
y(1:L-2, :) = y(1:L-2, :) + y(N+[1:L-2], :);
y = y(1:N, :);
y = cshift2D(y, 1-L/2);

if d == 2
   y = y';
end


function [u, v] = pm(a,b)

% [u v] = pm(a,b)
% u = (a + b)/sqrt(2);
% v = (a - b)/sqrt(2);

u = (a + b)/sqrt(2);
v = (a - b)/sqrt(2);

