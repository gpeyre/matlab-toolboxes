function y = perform_windowed_dct4_transform(x, w, dir, options)

% perform_windowed_dct4_transform - orthogonal local DCT
%
% y = perform_windowed_dct4_transform(x, w, dir, options);
%
%   This is an orthogonal transform, with x2 overlap, 
%   note that the low frequencies are removes if options.remove_lowfreq = 1.
%
%   Copyright (c) 2007 Gabriel Peyre


options.null = 0;
window_type = getoptions(options, 'window_type', 'Sine');
remove_lowfreq = getoptions(options, 'remove_lowfreq', 1);
    
if dir==1
    y = FastLDCT2ivAnalysis(x,window_type,w);
    y = y(2).coeff;
else
    c(1).coeff = zeros( size(x)/w );
    c(2).coeff = x;
    c(1).winwidth = w;
    c(2).winwidth = w;
    y = FastLDCT2ivSynthesis(c,window_type,w,remove_lowfreq);
end




function img = FastLDCT2ivSynthesis(coef,bellname,w,remove_low_freq) % pars3)

% FastLDCT2Synthesis -- 2D Local inverse DCT iv transform
%  Usage
%    img = FastLDCT2Synthesis(ldct,w)
%  Inputs
%    coef   	2D Local DCT structure array
%    w        	width of window
%    bellname name of bell to use, defaults to 'Sine'
%  Outputs
%    img	    2D reconstructed n by n image
%  Description
%    The matrix img contains image reconstructed from the Local DCT Decomposition.
% See Also
%   FastLDCT2Analysis, idct2
%

if nargin < 3 | bellname==0,
    bellname = 'Sine';
end

[n,J] = quadlength(coef(2).coeff);

if coef(2).winwidth ~= w
    error('Window width is different from given argument.');
    return;
end


img = zeros(n,n);
%
nbox = floor(n/w);
lfign=0.25;
if remove_low_freq
for boxcnt1=0:nbox-1
    for boxcnt2=0:nbox-1
        coef(2).coeff(boxcnt1*w+1:boxcnt1*w+1+floor(w*lfign),boxcnt2*w+1:boxcnt2*w+1+floor(w*lfign)) = 0;
    end
end
end

for ncol=1:n
    ldct = [struct('winwidth', w, 'coeff', 0) ...
        struct('winwidth', w, 'coeff', coef(2).coeff(:,ncol))];
    x = FastLDCTivSynthesis(ldct,bellname,w);
    img(:,ncol) = x;
end

for nrow=1:n
    ldct = [struct('winwidth', w, 'coeff', 0) ...
        struct('winwidth', w, 'coeff', img(nrow,:))];
    x = FastLDCTivSynthesis(ldct,bellname,w);
    img(nrow,:) = x';
end



%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%


function sig = FastLDCTivSynthesis(ldct,bellname,w,par3)
% FastLDCTivSynthesis -- Synthesize signal from local DCT iv coefficients (orthogonal fixed folding)
%  Usage
%    sig = FastLDCTivSynthesis(ldct,bellname,w)
%  Inputs
%    ldct       local DCT iv coefficients (structure array)
%    w		width of window
%    bell       name of bell to use, defaults to 'Sine'
%  Outputs
%    sig        signal whose orthonormal local DCT iv coeff's are ldct
%
%  See Also
%   FastLDCTivAnalysis, CPAnalysis, FCPSynthesis, fold, unfold, dct_iv, packet
%
[n,J] = dyadlength(ldct(2).coeff);
d = floor(log2(n/w));
%
% Create Bell
%
if nargin < 3,
    bellname = 'Sine';
end
m = n / 2^d /2;
[bp,bm] = MakeONBell(bellname,m);
%
%
%
x = zeros(1,n);
for b=0:(2^d-1),
    c = ldct(2).coeff(packet(d,b,n));
    y = dct_iv(c);
    [xc,xl,xr] = unfold(y,bp,bm);
    x(packet(d,b,n)) = x(packet(d,b,n)) + xc;
    if b>0,
        x(packet(d,b-1,n)) = x(packet(d,b-1,n)) + xl;
    else
        x(packet(d,0,n))   = x(packet(d,0,n)) + edgeunfold('left',xc,bp,bm);
    end
    if b < 2^d-1,
        x(packet(d,b+1,n)) = x(packet(d,b+1,n)) + xr;
    else
        x(packet(d,b,n))   = x(packet(d,b,n)) + edgeunfold('right',xc,bp,bm);
    end
end
sig = ShapeAsRow(x)';


%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%


function coef = FastLDCT2ivAnalysis(img,bellname,w,pars3)
% FastLDCT2ivAnalysis -- Analyze image into 2-d cosine packet coefficients at a given depth (window width)
%  Usage
%    coef = FastLDCT2ivAnalysis(img,w)
%  Inputs
%    img      2-d image to be transformed into basis
%    w        width of window
%    bellname name of bell to use, defaults to 'Sine'
%  Outputs
%    coef     2-d Local DCT iv coeffts
%
%  Description
%    Once a cosine packet basis has been selected (at a given depth),
%    this function may be used to expand a given
%    image in that basis.
%

if nargin < 3 | bellname==0,
    bellname = 'Sine';
end

[n,J] = quadlength(img);

d = floor(log2(n/w));

%
% CP image at depth d
%
coef = [struct('winwidth', w, 'coeff', zeros(floor(n/w),floor(n/w))) ...
    struct('winwidth', w, 'coeff', zeros(n,n))];
%
for nrow=1:n
    ldct = FastLDCTivAnalysis(img(nrow,:),bellname,w);
    coef(2).coeff(nrow,:) = ldct(2).coeff;
end

for ncol=1:n
    ldct = FastLDCTivAnalysis(coef(2).coeff(:,ncol),bellname,w);
    coef(2).coeff(:,ncol) = ldct(2).coeff;
end



function ldct = FastLDCTivAnalysis(x,bellname,w,par3)
% FastLDCTivAnalysis -- Local DCT iv transform (orthogonal fixed folding)
%  Usage
%    ldct = FastLDCTivAnalysis(x,bell,w)
%  Inputs
%    x        1-d signal:  length(x)=2^J
%    w        width of window
%    bell     name of bell to use, defaults to 'Sine'
%  Outputs
%    ldct     1-d Local DCT iv coefficients (structure array)
%  Description
%    The vector ldct contains coefficients of the Local DCT Decomposition.
% See Also
%   FastLDCTivSynthesis, CPAnalysis, FCPSynthesis, fold, unfold, dct_iv, packet
%

if nargin < 3,
    bellname = 'Sine';
end
[n,J] = dyadlength(x);

d = floor(log2(n/w));
%
% taper window
%
m = n / 2^d /2;
[bp,bm] = MakeONBell(bellname,m);
%
% packet table
%
n  = length(x);
ldct = [struct('winwidth', w, 'coeff', zeros(floor(n/w),1)) ...
    struct('winwidth', w, 'coeff', zeros(n,1))];
x  = ShapeAsRow(x);
%
nbox = 2^d;
for b=0:(nbox-1)
    if(b == 0) ,                         % gather packet and
        xc = x(packet(d,b,n));           % left, right neighbors
        xl = edgefold('left',xc,bp,bm);  % taking care of edge effects
    else
        xl = xc;
        xc = xr;
    end
    if (b+1 < nbox)
        xr = x(packet(d,b+1,n));
    else
        xr = edgefold('right',xc,bp,bm);
    end
    y = fold(xc,xl,xr,bp,bm);    % folding projection
    c = dct_iv(y);               % DCT-IV
    ldct(2).coeff(packet(d,b,n)) = c';  % store
end



%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%

function c = dct_iv(x)
% dct_iv -- Type (IV) Discrete Cosine Xform
%  Usage
%    c = dct_iv(x)
%  Inputs
%    x     1-d signal, length(x) = 2^J
%  Outputs
%    c     1-d cosine transform, length(x)=2^J
%
%  Description
%    The form c = dct_iv(x) computes c defined by
%         c_m = sqrt(2/N) * sum_n x(n) cos( pi * (m-.5) * (n-.5) / N )
%     where
%         1 <= m,n <= N,  N = length(x) = length(c)
%
%    To reconstruct, use the same function:
%         x = dct_iv(c)
%
%  See Also
%    CPAnalysis, CPSynthesis
%

n2 = 2*length(x);
y = zeros(1, 4*n2);
y(2:2:n2) = x(:);
z = fft(y);
c = sqrt(4/n2) .* real(z(2:2:n2));

%
% Copyright (c) 1993. David L. Donoho
%




%
%  Part of Wavelab Version 850
%  Built Tue Jan  3 13:20:40 EST 2006
%  This is Copyrighted Material
%  For Copying permissions see COPYING.m
%  Comments? e-mail wavelab@stat.stanford.edu

function [n,J] = dyadlength(x)
% dyadlength -- Find length and dyadic length of array
%  Usage
%    [n,J] = dyadlength(x)
%  Inputs
%    x    array of length n = 2^J (hopefully)
%  Outputs
%    n    length(x)
%    J    least power of two greater than n
%
%  Side Effects
%    A warning is issued if n is not a power of 2.
%
%  See Also
%    quadlength, dyad, dyad2ix
%
n = length(x) ;
J = ceil(log(n)/log(2));
if 2^J ~= n ,
    disp('Warning in dyadlength: n != 2^J')
end




%
%  Part of Wavelab Version 850
%  Built Tue Jan  3 13:20:40 EST 2006
%  This is Copyrighted Material
%  For Copying permissions see COPYING.m
%  Comments? e-mail wavelab@stat.stanford.edu


function [bp,bm] = MakeONBell(Name,m)
% MakeONBell -- Make Bell for Orthonormal Local Cosine Analysis
%  Usage
%    [bp,bm] = MakeONBell(bell,m)
%  Inputs
%    bell      bellname, currently 'Trivial','Sine'
%    m         length of bell
%  Outputs
%    bp        part of bell interior to domain
%    bm        part of bell exterior to domain
%
% See Also
%    CPAnalysis, CPSynthesis, MakeCosinePacket
%

Name = lower(Name);

xi = (1 + (.5:(m-.5))./m)./2;
if strcmp(Name,'trivial') || strcmp(Name,'constant')
    bp = sqrt(xi);
elseif strcmp(Name,'sine') || strcmp(Name,'sin')
    bp = sin( pi/2 .* xi );
else
    error('Unknown bell shape.');
end
bm = sqrt(1 - bp .^2);

%
% Copyright (c) 1993. David L. Donoho
%





%
%  Part of Wavelab Version 850
%  Built Tue Jan  3 13:20:40 EST 2006
%  This is Copyrighted Material
%  For Copying permissions see COPYING.m
%  Comments? e-mail wavelab@stat.stanford.edu


function row = ShapeAsRow(sig)
% ShapeAsRow -- Make signal a row vector
%  Usage
%    row = ShapeAsRow(sig)
%  Inputs
%    sig     a row or column vector
%  Outputs
%    row     a row vector
%
%  See Also
%    ShapeLike
%
row = sig(:)';





%
%  Part of Wavelab Version 850
%  Built Tue Jan  3 13:20:43 EST 2006
%  This is Copyrighted Material
%  For Copying permissions see COPYING.m
%  Comments? e-mail wavelab@stat.stanford.edu


function p = packet(d,b,n)
% packet -- Packet table indexing
%  Usage
%    p = packet(d,b,n)
%  Inputs
%    d     depth of splitting in packet decomposition
%    b     block index among 2^d possibilities at depth d
%    n     length of signal
%  Outputs
%    p     linear indices of all coeff's in that block
%

npack = 2^d;
p =  ( (b * (n/npack) + 1) : ((b+1)*n/npack ) ) ;

%
% Copyright (c) 1993. David L. Donoho
%




%
%  Part of Wavelab Version 850
%  Built Tue Jan  3 13:20:40 EST 2006
%  This is Copyrighted Material
%  For Copying permissions see COPYING.m
%  Comments? e-mail wavelab@stat.stanford.edu

function extra = edgefold(which,xc,bp,bm)
% edgefold -- Perform folding projection with (+,-) polarity at EDGES
%  Usage
%    extra = edgefold(which,xc,bp,bm)
%  Inputs
%    which  string, 'left'/'right', indicating which edge we are at
%    xc     unfolded data, center window
%    bp     interior of window
%    bm     exterior of window
%  Outputs
%    extra  pseudo-left/right packet
%
%  Description
%    The result should be used as either left or right packet in
%    fold to ensure exact reconstruction at edges.
%
%  See Also
%    fold, unfold, CPSynthesis, CPImpulse
%

n = length(xc);
m = length(bp);
back  = n:-1:(n-m+1);
front = 1:m;
extra = xc.*0;
%
if strcmp(which,'left'),
    extra(back) = xc(front) .* (1-bp)./bm;
else
    extra(front) = -xc(back) .* (1-bp)./bm;
end

%
%  Copyright (c) 1994. David L. Donoho
%




%
%  Part of Wavelab Version 850
%  Built Tue Jan  3 13:20:40 EST 2006
%  This is Copyrighted Material
%  For Copying permissions see COPYING.m
%  Comments? e-mail wavelab@stat.stanford.edu


function y = fold(xc,xl,xr,bp,bm)
% fold -- Folding projection with (+,-) polarity
%  Usage
%    y = fold(xc,xl,xr,bp,bm)
%  Inputs
%    xc,xl,xr    1-d signals: center, left, and right packets
%    bp,bm       interior and exterior of taper window
%  Outputs
%    y           the folded series, whose DCT gives the LCT coeffs
%
%  See Also
%    CPAnalysis
%
%  References
%    H. Malvar, IEEE ASSP, 1990.
%
%    Auscher, Weiss and Wickerhauser, in ``Wavelets: A Tutorial in
%    Theory and Applications,'' C. Chui, ed., Academic Press, 1992.
%

m = length(bp);
n = length(xc);
front = 1:m;
back  = n:-1:(n+1-m);
y = xc;
y(front)  = bp .* y(front)   + bm .* xl(back );
y(back )  = bp .* y(back )   - bm .* xr(front);

%
% Copyright (c) 1993. David L. Donoho
%




%
%  Part of Wavelab Version 850
%  Built Tue Jan  3 13:20:40 EST 2006
%  This is Copyrighted Material
%  For Copying permissions see COPYING.m
%  Comments? e-mail wavelab@stat.stanford.edu


function [xc,xl,xr] = unfold(y,bp,bm)
% unfold -- Undo folding projection with (+,-) polarity
%  Usage
%    [xc,xl,xr] = unfold(y,bp,bm)
%  Inputs
%    y     folded series
%    bp    interior of window
%    bm    exterior of window
%  Outputs
%    xc    unfolded central packet
%    xl    contribution of unfolding to left  packet
%    xr    contribution of unfolding to right packet
%
% See Also
%    fold, CPAnalysis, CPSynthesis, CPImpulse
%

	n = length(y);
	m = length(bp);
	xc = y;
	xl = 0 .*y;
	xr = 0 .*y;
	front = 1:m;
	back  = n:-1:(n+1-m);
	xc(front) =       bp .* y(front);
	xc(back)  =       bp .* y(back );
	xl(back)  =       bm .* y(front);
	xr(front) =     (-1) .* bm .* y(back);

%
% Copyright (c) 1993. David L. Donoho
%     
    
    
    
 
 
%
%  Part of Wavelab Version 850
%  Built Tue Jan  3 13:20:40 EST 2006
%  This is Copyrighted Material
%  For Copying permissions see COPYING.m
%  Comments? e-mail wavelab@stat.stanford.edu 
function extra = edgeunfold(which,xc,bp,bm)
% edgeunfold -- Undo folding projection with (+,-) polarity at EDGES
%  Usage
%    extra = endfold(which,xc,bp,bm)
%  Inputs
%    which  string, 'left'/'right', indicating which edge we are at
%    xc     unfolded data, center window, produced by unfold
%    bp     interior of window
%    bm     exterior of window
%  Outputs
%    extra  extra contribution to central packet
%
%  Description
%    The result should be added to the central packet to
%    ensure exact reconstruction at edges.
%
%  See Also
%    unfold, CPSynthesis, CPImpulse
%

	n = length(xc);
	m = length(bp);
	extra = xc.*0;
%
	if strcmp(which,'left'),
		front = 1:m;
		extra(front) = xc(front) .* (1-bp)./bp;
	else
		back  = n:-1:(n+1-m);
		extra(back) = xc(back) .* (1-bp)./bp;
	end	

%
%  Copyright (c) 1994. David L. Donoho
%
    
    
 
 
%
%  Part of Wavelab Version 850
%  Built Tue Jan  3 13:20:40 EST 2006
%  This is Copyrighted Material
%  For Copying permissions see COPYING.m
%  Comments? e-mail wavelab@stat.stanford.edu 



function [n,J] = quadlength(x)
% quadlength -- Find length and dyadic length of square matrix
%  Usage
%    [n,J] = quadlength(x)
%  Inputs
%    x   2-d image; size(n,n), n = 2^J (hopefully)
%  Outputs
%    n   length(x)
%    J   least power of two greater than n
%
%  Side Effects
%    A warning message is issue if n is not a power of 2,
%    or if x is not a square matrix.
%
s = size(x);
n = s(1);
if s(2) ~= s(1),
    disp('Warning in quadlength: nr != nc')
end
k = 1 ; J = 0; while k < n , k=2*k; J = 1+J ; end ;
if k ~= n ,
    disp('Warning in quadlength: n != 2^J')
end

%
% Copyright (c) 1993. David L. Donoho
%
