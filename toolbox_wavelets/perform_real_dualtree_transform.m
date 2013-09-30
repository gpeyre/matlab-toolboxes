function y = perform_real_dualtree_transform(x,Jmin,options)

% perform_real_dualtree_transform - perform a real valued dual-tree wavelet transform
%
%   y = perform_real_dualtree_transform(x,Jmin,options);
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
    y1 = dualtree2D(x, J, Faf, af);
    % flatten everything
    y = {};
    for j=1:J
        for d1=1:2
            for d2=1:3
                y{end+1} = y1{j}{d1}{d2};
            end
        end
    end
    for d1=1:2
        y{end+1} = y1{end}{d1};
    end
else
    % flatten everything
    x1 = {};
    k = 0;
    for j=1:J
        for d1=1:2
            for d2=1:3
                k = k+1;
                x1{j}{d1}{d2} = x{k};
            end
        end
    end
    for d1=1:2
        k = k+1;
        x1{J+1}{d1} = x{k};
    end
    y = idualtree2D(x1, J, Fsf, sf);
end



function w = dualtree2D(x, J, Faf, af)

% 2D Dual-Tree Discrete Wavelet Transform
%
% USAGE:
%   w = dualtree2D(x, J, Faf, af)
% INPUT:
%   x - M by N array
%   J - number of stages
%   Faf - first stage filters
%   af - filters for remaining stages
% OUPUT:
%   w{j}{d1}{d2} - DWT coefficients
%       j = 1..J, k = 1..2, d = 1..3
%   w{J+1}{k} - lowpass coefficients
%       k = 1..2
% % EXAMPLE:
%   x = rand(256);
%   J = 3;
%   [Faf, Fsf] = FSfarras;
%   [af, sf] = dualfilt1;
%   w = dualtree2D(x, J, Faf, af);
%   y = idualtree2D(w, J, Fsf, sf);
%   err = x - y;
%   max(max(abs(err)))
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

% normalization
x = x/sqrt(2);

% Tree 1
[x1 w{1}{1}] = afb2D(x, Faf{1});      % stage 1
for j = 2:J
    [x1 w{j}{1}] = afb2D(x1, af{1});  % remaining stages
end
w{J+1}{1} = x1;                       % lowpass subband

% Tree 2
[x2 w{1}{2}] = afb2D(x, Faf{2});      % stage 1
for j = 2:J
    [x2 w{j}{2}] = afb2D(x2, af{2});  % remaining stages
end
w{J+1}{2} = x2;                       % lowpass subband

% sum and difference
for j = 1:J
    for m = 1:3
        A = w{j}{1}{m};
        B = w{j}{2}{m};
        w{j}{1}{m} = (A+B)/sqrt(2);
        w{j}{2}{m} = (A-B)/sqrt(2);
    end
end



function y = idualtree2D(w, J, Fsf, sf)

% Inverse 2-D Dual-Tree Discrete Wavelet Transform
% 
% USAGE:
%   y = idualtree2D(w, J, Fsf, sf)
% INPUT:
%   J - number of stages
%   Fsf - synthesis filters for final stage
%   sf -  synthesis filters for preceeding stages
% OUPUT:
%   y - output array
% See idualtree2D
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

% sum and difference
for k = 1:J
    for m = 1:3
        A = w{k}{1}{m};
        B = w{k}{2}{m};
        w{k}{1}{m} = (A+B)/sqrt(2);
        w{k}{2}{m} = (A-B)/sqrt(2);
    end
end

% Tree 1
y1 = w{J+1}{1};
for j = J:-1:2
   y1 = sfb2D(y1, w{j}{1}, sf{1});
end
y1 = sfb2D(y1, w{1}{1}, Fsf{1});

% Tree 2
y2 = w{J+1}{2};
for j = J:-1:2
   y2 = sfb2D(y2, w{j}{2}, sf{2});
end
y2 = sfb2D(y2, w{1}{2}, Fsf{2});

% normalization
y = (y1 + y2)/sqrt(2);



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

