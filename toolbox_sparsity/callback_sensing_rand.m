function y = callback_sensing_rand(x, dir, options)

% callback_sensing_rand - perform random sensing
%
%   y = callback_sensing_rand(x, dir, options);
%
%   compute y=K*x (dir=1) or y=K^{*}*x (dir=-1) or y=K^{+}*x (pseudo inverse) 
%   where K is a random matrix. 
%
%   You need to set options.n and options.p the dimension
%   of the (p,n) matrix K (P=#measurements)
%
%   options.cs_type selects different kinds of matrix
%       'user' : you provide your own matrix
%       'fourier' : random sub-selection of frequencies
%       'fourier-scramb' : random sub-selection of frequencies of scrambled signal
%       'hadamard' : random sub-selection of hadamard projections
%       'fourier-scramb' : random sub-selection of hadamard projections of scrambled signal
%       'noiselets' : random sub-selection of noiselets projections
%       'noiselets-scramb' : random sub-selection of noiselets projections of scrambled signal
%       'sinus' : random sub-selection of discrete sinus transform projections
%       'sinus-scramb' : random sub-selection of discrete sinus transform projectionsof scrambled signal
%
% options.rand_matrix (for  using your own matrix)
% or options.n (size of the signal) and options.p 
% (number of measurements) for automatic generation.
%
%   Copyright (c) 2008 Gabriel Peyre


options.null = 0;
cs_type = getoptions(options, 'cs_type', 'sinus');
cs_matrix = getoptions(options, 'cs_matrix', []);

csndims = getoptions(options, 'csndims', 1);

if csndims==2
    options.csndims = 1;
    if dir==1
        y = callback_sensing_rand(x(:), dir, options);
    else
        y = callback_sensing_rand(x(:), dir, options);
        y = reshape(y, sqrt(length(y))*[1 1]);
    end
    return;
end

if not(isempty(cs_matrix))
    n = size(cs_matrix, 2);
    p = size(cs_matrix, 1);
else
    n = getoptions(options, 'n', 0, 1);
    p = getoptions(options, 'p', 0, 1);
end

switch lower(cs_type)
    case {'sinus' 'hadamard' 'fourier' 'sinus-scramb' 'hadamard-scramb' 'fourier-scramb' }
        y = cs_operator(x, p,n, dir, cs_type);


    case 'user'
        if isempty(cs_matrix)
            error('You must provide options.cs_matrix');
        end
        
        if dir==1
            y = cs_matrix*x;
        elseif dir==-1
            y = cs_matrix'*x;
        elseif dir==-2
            cs_matrix_pinv = getoptions(options, 'cs_matrix_pinv', []);
            if isempty(cs_matrix_pinv)
                warning('You should provide options.cs_matrix_pinv');
                cs_matrix_pinv = pinv(cs_matrix);
            end
            y = cs_matrix_pinv*x;
        else
            error('dir should be +1,-1 or -2');
        end

    otherwise 
        error('Unknown CS operator');
end



function y = cs_operator(x, m,n, dir, type)

% FastCSOperator: The operator form of a random sampling matrix for the 
% compressed sensing problem.
% Specifically, it returns y = A(:,I)*x (mode = 1) or y = A(:,I)'*x (mode = 2),
% where A is an mxdim random sampling matrix defined as
% A = P*H*Q, where P,Q are random permutation matrices, 
% H is a fast Hadamard/Fourier operator, and I is 
% a subset of the columns of A, i.e. a subset of 1:dim of length n.

% Pstate - state of random generator for P matrix
% Qstate - state of random generator for Q matrix
Pstate = 4972169;
Qstate = 7256157;

do_scramb = 0;
if strcmp(type, 'sinus-scramb') || strcmp(type, 'hadamard-scramb') || strcmp(type, 'fourier-scramb')
    do_scramb = 1;
end

ntest = size(x,2);

if (dir == +1) % analysis
    % Apply matrix Q
    rand('state', Qstate);
    if do_scramb
        x = x( randperm(n),: );
    end
    % Apply matrix H
    switch type
        case {'sinus' 'sinus-scramb'}
            x = RST(x);
        case {'hadamard' 'hadamard-scramb'}
            x = FHT(x);
        case {'fourier' 'fourier-scramb'}
            x = fft(x)/sqrt(n);
        otherwise 
            error('Unknown CS operator');
    end
    % Apply matrix P
    rand('state', Pstate);
    sel = randperm(n);
    y = x(sel(1:m),:);
else % Adjoint operator
    % Apply matrix P^T
    rand('state', Pstate);
    sel = randperm(n);
    y = zeros(n,ntest);
    y( sel(1:m),: ) = x;
    % Apply matrix H^T
    switch type
        case {'sinus' 'sinus-scramb'}
            y = Inv_RST(y);
        case {'hadamard' 'hadamard'}
            y = Inv_FHT(y);
        case {'fourier' 'fourier-scramb'}
            y = ifft(y)*sqrt(n);
        otherwise 
            error('Unknown CS operator');
    end
    % Apply matrix Q^T
    if do_scramb
        rand('state', Qstate);
        y(randperm(n),:) = y;
    end
end


function S = RST(X)

% RST: Real Sinusoid Transform of an n-vector
%  Usage:
%    S = RST(X);
%  Inputs:
%    X      input n-vector 
%  Outputs:
%    S      output vector which contains the transform coeffs. 
%
% Description
%   RST computes a 1-D real sinusoid transform of an n vector, n dyadic,
%   by taking the fft of X and using conjugate symmetry to eliminate
%   complex values. 
% See Also
%   Inv_RST, RST2

n = size(X,1);
S = fft_mid0(X) ./ sqrt(n);

n2 = n/2 + 1;    % Center point

S(2:(n2-1),:) = sqrt(2) .* real(S(2:(n2-1),:));
S((n2+1):n,:) = -sqrt(2) .* imag(S((n2+1):n,:));

%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%

function X = Inv_RST(S)
% Inv_RST: Inverse Real Sinusoid Transform of an n vector
%  Usage:
%    X = Inv_RST(S);
%  Inputs:
%    S      n vector 
%  Outputs:
%    X      reconstructed vector. 
%
% Description
%   Inv_RST computes the inverse of RST by rearranging frequencies and 
%   taking the inverse 1-D fft.
% See Also
%   RST

n = size(S,1);
n2 = n/2 + 1;    % Center point
isqrt2 = 1./sqrt(2);

X = S;

% Top row
X(2:(n2-1),:)   = isqrt2 .* (S(2:(n2-1),:) + i .* S(n:-1:(n2+1),:));
X(n:-1:(n2+1),:) = isqrt2 .* (S(2:(n2-1),:) - i .* S(n:-1:(n2+1),:));

X = real(ifft_mid0(X)) .* sqrt(n);

%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%


function y=fft_mid0(x)
% fftmid0 -- 1d fft with argument [-pi,pi] (instead of [0,2pi])
% Usage:
%   Y = fft_mid0(X)
% Inputs:
%   X	Array(n) 
% Outputs:
%   Y   Array(n)
% Description:
%  Performs 1d fft with grid (-n/2):(n/2-1) on both time
%  and frequency side. 
%    y(k) = sum_{t=-n/2}^{n/2-1} exp(i 2pi/n kt) x(t) , (-n/2) <= k < n/2
%

y = fftshift(fft(fftshift(x,1)),1);

%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%


function y=ifft_mid0(x)
% ifft_mid0 -- 1d fft with argument [-pi,pi] (instead of [0,2pi])
% Usage:
%   Y = ifft_mid0(X)
% Inputs:
%   X	Array(n) 
% Outputs:
%   Y   Array(n)
% Description:
%  Performs 1d fft with grid (-n/2):(n/2-1) on both time
%  and frequency side. 
%    y(k) = sum_{t=-n/2}^{n/2-1} exp(i 2pi/n kt) x(t) , (-n/2) <= k < n/2
%

y = fftshift(ifft(fftshift(x,1)),1);

%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%

function y = FHT(x)
% FHT: Computes the Hadamard transform of a signal
%  Usage:
%    y = FHT(x);
%  Inputs:
%    x      signal of dyadic length
%  Outputs:
%    y      Hadamard transform coefficients
%
% Description
%   FHT computes the Fast Hadamard transform of the signal x.
% See Also
%   IFHT

%x = x(:);
n = size(x,1);

y = x;
t = zeros(n,size(x,2));  

odds = 1:2:n;
evens = 2:2:n;

for ii = 1:log2(n)
    t(odds,:) = y(odds,:) + y(evens,:);
    t(evens,:) = y(odds,:) - y(evens,:);
    y(1:(n/2),:) = t(odds,:);
    y((n/2+1):n,:) = t(evens,:);
end

y = y ./ sqrt(n);
%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%

function y = Inv_FHT(x)
% Inv_FHT: Computes the inverse Hadamard transform of a signal
%  Usage:
%    y = Inv_FHT(x);
%  Inputs:
%    x      Hadamard transform coefficients
%  Outputs:
%    y      signal
%
% Description
%   Inv_FHT applies the forward Hadamard transform, since it is a self-adjoint operator.
% See Also
%   FHT

y = FHT(x);
%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%

