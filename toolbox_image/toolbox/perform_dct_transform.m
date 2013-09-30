function y = perform_dct_transform(x,dir)

% perform_dct_transform - discrete cosine transform
%
%   y = perform_dct_transform(x,dir);
%
%   Copyright (c) 2006 Gabriel Peyre



if size(x,1)==1 || size(x,2)==1
    % 1D transform
    if dir==1
        y = dct(x);
    else
        y = idct(x);
    end
else
    if dir==1
        y = dct2(x);
        %   y = perform_dct2_transform(x);
    else
        y = idct2(x);
    end
end




function b=dct(a,n)
%DCT  Discrete cosine transform.
%
%   Y = DCT(X) returns the discrete cosine transform of X.
%   The vector Y is the same size as X and contains the
%   discrete cosine transform coefficients.
%
%   Y = DCT(X,N) pads or truncates the vector X to length N
%   before transforming.
%
%   If X is a matrix, the DCT operation is applied to each
%   column.  This transform can be inverted using IDCT.
%
%   See also FFT, IFFT, IDCT.

%   Author(s): C. Thompson, 2-12-93
%              S. Eddins, 10-26-94, revised
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.7 $  $Date: 2002/04/15 01:10:40 $

%   References:
%   1) A. K. Jain, "Fundamentals of Digital Image
%      Processing", pp. 150-153.
%   2) Wallace, "The JPEG Still Picture Compression Standard",
%      Communications of the ACM, April 1991.


if nargin == 0,
    error('Not enough input arguments.');
end

if isempty(a)
    b = [];
    return
end

% If input is a vector, make it a column:
do_trans = (size(a,1) == 1);
if do_trans, a = a(:); end

if nargin==1,
    n = size(a,1);
end
m = size(a,2);

% Pad or truncate input if necessary
if size(a,1)<n,
    aa = zeros(n,m);
    aa(1:size(a,1),:) = a;
else
    aa = a(1:n,:);
end

% Compute weights to multiply DFT coefficients
ww = (exp(-i*(0:n-1)*pi/(2*n))/sqrt(2*n)).';
ww(1) = ww(1) / sqrt(2);

if rem(n,2)==1 | ~isreal(a), % odd case
    % Form intermediate even-symmetric matrix
    y = zeros(2*n,m);
    y(1:n,:) = aa;
    y(n+1:2*n,:) = flipud(aa);

    % Compute the FFT and keep the appropriate portion:
    yy = fft(y);
    yy = yy(1:n,:);

else % even case
    % Re-order the elements of the columns of x
    y = [ aa(1:2:n,:); aa(n:-2:2,:) ];
    yy = fft(y);
    ww = 2*ww;  % Double the weights for even-length case
end

% Multiply FFT by weights:
b = ww(:,ones(1,m)) .* yy;

if isreal(a), b = real(b); end
if do_trans, b = b.'; end



function b=dct2(arg1,mrows,ncols)
%DCT2 Compute 2-D discrete cosine transform.
%   B = DCT2(A) returns the discrete cosine transform of A.
%   The matrix B is the same size as A and contains the
%   discrete cosine transform coefficients.
%
%   B = DCT2(A,[M N]) or B = DCT2(A,M,N) pads the matrix A with
%   zeros to size M-by-N before transforming. If M or N is
%   smaller than the corresponding dimension of A, DCT2 truncates
%   A.
%
%   This transform can be inverted using IDCT2.
%
%   Class Support
%   -------------
%   A can be numeric or logical. The returned matrix B is of
%   class double.
%
%   Example
%   -------
%       RGB = imread('autumn.tif');
%       I = rgb2gray(RGB);
%       J = dct2(I);
%       imshow(log(abs(J)),[]), colormap(jet), colorbar
%
%   The commands below set values less than magnitude 10 in the
%   DCT matrix to zero, then reconstruct the image using the
%   inverse DCT function IDCT2.
%
%       J(abs(J)<10) = 0;
%       K = idct2(J);
%       imview(I)
%       imview(K,[0 255])
%
%   See also FFT2, IDCT2, IFFT2.

%   Copyright 1993-2003 The MathWorks, Inc.
%   $Revision: 5.22.4.2 $  $Date: 2003/05/03 17:50:23 $

%   References:
%        1) A. K. Jain, "Fundamentals of Digital Image
%           Processing", pp. 150-153.
%        2) Wallace, "The JPEG Still Picture Compression Standard",
%           Communications of the ACM, April 1991.

[m, n] = size(arg1);
% Basic algorithm.
if (nargin == 1),
    if (m > 1) & (n > 1),
        b = dct(dct(arg1).').';
        return;
    else
        mrows = m;
        ncols = n;
    end
end

% Padding for vector input.
a = arg1;
if nargin==2, ncols = mrows(2); mrows = mrows(1); end
mpad = mrows; npad = ncols;
if m == 1 & mpad > m, a(2, 1) = 0; m = 2; end
if n == 1 & npad > n, a(1, 2) = 0; n = 2; end
if m == 1, mpad = npad; npad = 1; end   % For row vector.

% Transform.

b = dct(a, mpad);
if m > 1 & n > 1, b = dct(b.', npad).'; end



function a = idct(b,n)
%IDCT Inverse discrete cosine transform.
%
%   X = IDCT(Y) inverts the DCT transform, returning the original
%   vector if Y was obtained using Y = DCT(X).
%
%   X = IDCT(Y,N) pads or truncates the vector Y to length N
%   before transforming.
%
%   If Y is a matrix, the IDCT operation is applied to each
%   column.
%
%   See also FFT,IFFT,DCT.

%   Copyright 1993-2003 The MathWorks, Inc.
%   $Revision: 5.12.4.1 $  $Date: 2003/01/26 05:59:37 $

%   References:
%       1) A. K. Jain, "Fundamentals of Digital Image
%          Processing", pp. 150-153.
%       2) Wallace, "The JPEG Still Picture Compression Standard",
%          Communications of the ACM, April 1991.

% checknargin(1,2,nargin,mfilename);

if ~isa(b, 'double')
    b = double(b);
end

if min(size(b))==1
    if size(b,2)>1
        do_trans = 1;
    else
        do_trans = 0;
    end
    b = b(:);
else
    do_trans = 0;
end
if nargin==1,
    n = size(b,1);
end
m = size(b,2);

% Pad or truncate b if necessary
if size(b,1)<n,
    bb = zeros(n,m);
    bb(1:size(b,1),:) = b;
else
    bb = b(1:n,:);
end

if rem(n,2)==1 | ~isreal(b), % odd case
    % Form intermediate even-symmetric matrix.
    ww = sqrt(2*n) * exp(j*(0:n-1)*pi/(2*n)).';
    ww(1) = ww(1) * sqrt(2);
    W = ww(:,ones(1,m));
    yy = zeros(2*n,m);
    yy(1:n,:) = W.*bb;
    yy(n+2:n+n,:) = -j*W(2:n,:).*flipud(bb(2:n,:));

    y = ifft(yy);

    % Extract inverse DCT
    a = y(1:n,:);

else % even case
    % Compute precorrection factor
    ww = sqrt(2*n) * exp(j*pi*(0:n-1)/(2*n)).';
    ww(1) = ww(1)/sqrt(2);
    W = ww(:,ones(1,m));

    % Compute x tilde using equation (5.93) in Jain
    y = ifft(W.*bb);

    % Re-order elements of each column according to equations (5.93) and
    % (5.94) in Jain
    a = zeros(n,m);
    a(1:2:n,:) = y(1:n/2,:);
    a(2:2:n,:) = y(n:-1:n/2+1,:);
end

if isreal(b), a = real(a); end
if do_trans, a = a.'; end

function a = idct2(arg1,mrows,ncols)
%IDCT2 Compute 2-D inverse discrete cosine transform.
%   B = IDCT2(A) returns the two-dimensional inverse discrete
%   cosine transform of A.
%
%   B = IDCT2(A,[M N]) or B = IDCT2(A,M,N) pads A with zeros (or
%   truncates A) to create a matrix of size M-by-N before
%   transforming.
%
%   For any A, IDCT2(DCT2(A)) equals A to within roundoff error.
%
%   The discrete cosine transform is often used for image
%   compression applications.
%
%   Class Support
%   -------------
%   The input matrix A can be of class double or of any
%   numeric class. The output matrix B is of class double.
%
%   See also DCT2, DCTMTX, FFT2, IFFT2.

%   Copyright 1993-2003 The MathWorks, Inc.
%   $Revision: 5.17.4.1 $  $Date: 2003/01/26 05:55:39 $

%   References:
%   1) A. K. Jain, "Fundamentals of Digital Image
%      Processing", pp. 150-153.
%   2) Wallace, "The JPEG Still Picture Compression Standard",
%      Communications of the ACM, April 1991.

% checknargin(1,3,nargin,mfilename);

[m, n] = size(arg1);
% Basic algorithm.
if (nargin == 1),
    if (m > 1) & (n > 1),
        a = idct(idct(arg1).').';
        return;
    else
        mrows = m;
        ncols = n;
    end
end

% Padding for vector input.

b = arg1;
if nargin==2,
    ncols = mrows(2);
    mrows = mrows(1);
end

mpad = mrows; npad = ncols;
if m == 1 & mpad > m, b(2, 1) = 0; m = 2; end
if n == 1 & npad > n, b(1, 2) = 0; n = 2; end
if m == 1, mpad = npad; npad = 1; end   % For row vector.

% Transform.

a = idct(b, mpad);
if m > 1 & n > 1, a = idct(a.', npad).'; end
