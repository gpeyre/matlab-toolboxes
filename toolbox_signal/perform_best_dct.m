function [y,btree] = perform_best_dct(x,dir,options)

% perform_best_dct - compute a best local DCT.
%
%   [y,btree] = perform_best_dct(x,dir,options);
%
%   This code is based on the wavelab implementation.
%   You can give the basis parameter in options.btree.
%
%   Copyright (c) 2006 Gabriel Peyré

global WLVERBOSE;
WLVERBOSE = 'False';

bell = getoptions(options, 'bell', 'Sine');

options.null = 0;
par = [];
if isfield(options, 'ent')
    ent = options.ent;
else
    ent = 'Entropy';
    ent = 'lagrangian';
    ent = 'l^p'; par = 1;
    ent = 'l1';
end
mu0 = getoptions(options, 'mu0', 7);
Jmin = getoptions(options, 'Jmin', 1);
btree = getoptions(options, 'btree', []);


fixedwidth = getoptions(options, 'fixedwidth', -1);

if fixedwidth>0
    if dir==1
        y = FastLDCTivAnalysis(x,bell,fixedwidth);
    else
        y = FastLDCTivSynthesis(x,bell,fixedwidth);
    end
    btree = [];
    return;
end

x = x(:);
n = size(x,1);
Jmax = log2(n)-1;
J = Jmax-Jmin+1;

if ~isempty(btree) && btree(1)==-1
    % special meaning for DCT implementation
    if dir==1
        y = dct(x);
    else
        y = idct(x);
    end
    return;
end

if isempty(btree)
    if isfield(options, 'T')
        T = options.T;
    else
        T = 10;
    end
    % compute best tree
    if dir~=1
        error('You must provide btree for backward transform.');
    end
    % perform local DCT
    cp = CPAnalysis(x,J,bell);
    % compute local lagrangian
    switch lower(ent)
        case 'lagrangian'
            [e0,e1,e2] = compute_lagrangian_dct(cp,T);
            % compute penalization
            pen = (e0*0 + 1)/mu0;
            stree = e2 + T^2*e0 + T^2*pen;
        case 'l1'
            [e0,e1,e2] = compute_lagrangian_dct(cp,T);
            % compute penalization
            pen = (e0*0 + 1)/mu0;
            stree = T*e1 + T^2*pen;
        otherwise
            stree = CalcStatTree(cp,ent,par) * norm(cp(:,1));
    end
    % compute best tree
    [btree,vtree] = BestBasis(stree,J);
    if 0
        PlotPacketTable(cp);
        PlotBasisTree(btree,J,stree,'');
    end
end

if dir==1
    %%% forward %%%
    y = fpt_cp(btree,x,J,bell);
    % ImagePhasePlane('CP',btree,cp,'');
else
    %%% forward transform
    y = ipt_cp(btree,x,J,bell);
end
y = y(:);


%%% non adaptative

function coef = FastLDCTivAnalysis(x,bellname,w,par3)
% FastLDCTivAnalysis -- Local DCT iv transform (orthogonal fixed folding)
%  Usage
%    ldct = FastLDCTivAnalysis(x,bell,w)
%  Inputs
%    x        1-d signal:  length(x)=2^J
%    w        width of window
%    bell     name of bell to use, defaults to 'Sine'
%  Outputs
%    coef     1-d Local DCT iv coefficients
%  Description
%    The vector coef contains coefficients of the Local DCT Decomposition.
% See Also
%   FastLDCTivSynthesis, CPAnalysis, FCPSynthesis, fold, unfold, dct_iv, packet
%

if nargin < 3 | bellname==0,
    bellname = 'Sine';
end

[n,J] = dyadlength(x);

d = floor(log2(n/w));

%
% CP image at depth d
%
%

%
% taper window
%
m = n / 2^d /2;
[bp,bm] = MakeONBell(bellname,m);
%
% packet table
%
n  = length(x);
x  = ShapeAsRow(x);
coef = zeros(n,1);
%
nbox = 2^d;
for b=0:(nbox-1)
    if(b == 0) ,                             % gather packet and
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
    coef(packet(d,b,n)) = c';  % store
end



%%%

function x = FastLDCTivSynthesis(coef,bellname,w,pars3)
% FastLDCTivSynthesis -- Synthesize signal from local DCT iv coefficients (orthogonal fixed folding)
%  Usage
%    sig = FastLDCTivSynthesis(ldct,bellname,w)
%  Inputs
%    coef       local DCT iv coefficients
%    w		width of window
%    bell       name of bell to use, defaults to 'Sine'
%  Outputs
%    x          signal whose orthonormal local DCT iv coeff's are ldct
%
%  See Also
%   FastLDCTivAnalysis, CPAnalysis, FCPSynthesis, fold, unfold, dct_iv, packet
%

[n,J] = dyadlength(coef);
d = floor(log2(n/w));

%
% Create Bell
%
if nargin < 4 | bellname==0,
    bellname = 'Sine';
end
m = n / 2^d /2;
[bp,bm] = MakeONBell(bellname,m);

nbox = floor(n/w);
%lfign=0.25;
lfign=-1;
for boxcnt=0:nbox-1
    coef(boxcnt*w+1:boxcnt*w+1+floor(w*lfign)) = 0;
end
%
%
%
x = zeros(1,n);
for b=0:(2^d-1),
    c = coef(packet(d,b,n));
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

x = x(:);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [e0,e1,e2] = compute_lagrangian_dct(cp,T)

% compute_lagrangian_dct - compute norms for local basis
%
%   [e0,e1,e2] = compute_lagrangian_dct(cp,T);
%
% if v is the vector of dot product for the ith local DCT decomposition,
% then :
%   e0(i) = #{ i \ |v[i]|>T }
%   e1(i) = \sum |v[i]|
%   e2(i) = \sum_{ i \ |v[i]|<T } v[i]^2
%
%   Copyright (c) 2006 Gabriel Peyré

n = size(cp,1);
J = size(cp,2);
e0 = zeros(2^J-1,1);
e1 = zeros(2^J-1,1);
e2 = zeros(2^J-1,1);
s = 0;
for j=0:J-1
    for k=1:2^j
        s = s+1;
        sel = (k-1)*n/2^j+1:k*n/2^j;
        x = cp( sel ,j+1);
        e0(s) = sum( abs(x)>T );
        e1(s) = sum( abs(x) );
        I = find( abs(x)<T ); %  x(I) = 0;
        e2(s) = sum( x(I).^2 );
    end
end





function coef = fpt_cp(basis,x,D,bell)
% FPT_CP -- Fast transform into specific cosine packet basis
%  Usage
%    coef = FPT_CP(basis,x,D,bell)
%  Inputs
%    basis    tree selecting cosine packet basis
%    x        1-d signal to be transformed into basis
%    D,bell	  maximum depth of tree, type of bell
%  Outputs
%    coef     1-d cosine packet coeffts in given basis
%
%  Description
%    Once a cosine packet basis has been selected (presumably via
%    BestBasis), this function may be used to expand a given signal
%    in that basis.
%
[n,J] = dyadlength(x);
if nargin < 4,
    bell = 'Sine';
end
if nargin < 3,
    D = min(7,J-3);
end
coef = ShapeAsRow(x);


%   setup bell
m = n / 2^D /2;
[bp,bm] = MakeONBell(bell,m);

% At scale 0,  should fold around edges
% Dangling; to be added in a later version

% initialize tree traversal stack
stack = zeros(2,100); % column = [d, b]'

% pushdown root
k = 1;
stack(:,k) = [0 0]'; % d, b

while(k > 0),

    %  pop stack
    d = stack(1,k);   % depth of this node
    b = stack(2,k);   % branch at this depth
    k = k-1;
    %fprintf('d b'); disp([d b])

    if(basis(node(d,b)) ~= 0) ,

        % nonterminal node; fold around middle

        lo = 1 + b*n/2^d; hi = (b+1)*n/2^d;
        %fprintf('[lo hi]'); disp([lo hi])

        % fold middle
        midpost = floor((lo+hi)/2) + (1:m);
        midpre  = ceil ((lo+hi)/2) - (1:m);
        cf_right = coef(midpost);
        cf_left  = coef(midpre );
        coef(midpost) = bp .* cf_right + bm .* cf_left ;
        coef(midpre ) = bp .* cf_left  - bm .* cf_right;

        % pushdown children
        k = k+1; stack(:,k) = [(d+1) (2*b)   ]';
        k = k+1; stack(:,k) = [(d+1) (2*b+1) ]';

    else

        % terminal node -- analyze by dct_iv
        sig = coef(packet(d,b,n));
        coef(packet(d,b,n)) = dct_iv(sig);;

    end
end

%
% Copyright (c) 1993. David L. Donoho
%


%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%

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


function coef = FPT_CP(basis,x,D,bell)
% FPT_CP -- Fast transform into specific cosine packet basis
%  Usage
%    coef = FPT_CP(basis,x,D,bell)
%  Inputs
%    basis    tree selecting cosine packet basis
%    x        1-d signal to be transformed into basis
%    D,bell	  maximum depth of tree, type of bell
%  Outputs
%    coef     1-d cosine packet coeffts in given basis
%
%  Description
%    Once a cosine packet basis has been selected (presumably via
%    BestBasis), this function may be used to expand a given signal
%    in that basis.
%
[n,J] = dyadlength(x);
if nargin < 4,
    bell = 'Sine';
end
if nargin < 3,
    D = min(7,J-3);
end
coef = ShapeAsRow(x);


%   setup bell
m = n / 2^D /2;
[bp,bm] = MakeONBell(bell,m);

% At scale 0,  should fold around edges
% Dangling; to be added in a later version

% initialize tree traversal stack
stack = zeros(2,100); % column = [d, b]'

% pushdown root
k = 1;
stack(:,k) = [0 0]'; % d, b

while(k > 0),

    %  pop stack
    d = stack(1,k);   % depth of this node
    b = stack(2,k);   % branch at this depth
    k = k-1;
    %fprintf('d b'); disp([d b])

    if(basis(node(d,b)) ~= 0) ,

        % nonterminal node; fold around middle

        lo = 1 + b*n/2^d; hi = (b+1)*n/2^d;
        %fprintf('[lo hi]'); disp([lo hi])

        % fold middle
        midpost = floor((lo+hi)/2) + (1:m);
        midpre  = ceil ((lo+hi)/2) - (1:m);
        cf_right = coef(midpost);
        cf_left  = coef(midpre );
        coef(midpost) = bp .* cf_right + bm .* cf_left ;
        coef(midpre ) = bp .* cf_left  - bm .* cf_right;

        % pushdown children
        k = k+1; stack(:,k) = [(d+1) (2*b)   ]';
        k = k+1; stack(:,k) = [(d+1) (2*b+1) ]';

    else

        % terminal node -- analyze by dct_iv
        sig = coef(packet(d,b,n));
        coef(packet(d,b,n)) = dct_iv(sig);;

    end
end

%
% Copyright (c) 1993. David L. Donoho
%


%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%



function [basis,value] = BestBasis(tree,D)
% BestBasis -- Coifman-Wickerhauser Best-Basis Algorithm
%  Usage
%    [btree,vtree] = BestBasis(stree,D)
%  Inputs
%    stree    stat-tree (output by CalcStatTree)
%    D        maximum depth of tree-search
%  Outputs
%    btree    basis-tree of best basis
%    vtree    value of components of best basis
%             vtree(1) holds value of best basis
%
%  Description
%    The best-basis algorithm is used to pick out the ``best''
%    basis from all the possible bases in the packet table.
%    Here ``best'' means minimizing an additive measure of
%    information, called entropy by Coifman and Wickerhauser.
%
%    Once the stattree of entropy values is created, BestBasis
%    selects the best basis using the pruning algorithm described in
%    Wickerhauser's book.
%
%  Examples
%    [n,D] = dyadlength(signal);
%    qmf = MakeONFilter('Coiflet',3);
%    wp = WPAnalysis(signal,D,qmf);
%    stree = CalcStatTree(wp,'Entropy');
%    [btree,vtree] = BestBasis(stree,D);
%
%  Algorithm
%    Yale University has filed a patent application for this algorithm.
%    Commercial Development based on this algorithm should be cleared
%    by Yale University. Contact them for licensing information.
%
%  See Also
%    WPAnalysis, CalcStatTree, CPTour, WPTour
%
%  References
%    Wickerhauser, M.V. _Adapted_Wavelet_Analysis_
%
global WLVERBOSE
basis = zeros(size(tree));
value = tree;
for d=D-1:-1:0,
    for b=0:(2^d-1),
        vparent = tree(node(d,b));
        vchild  = value(node(d+1,2*b)) + value(node(d+1,2*b+1));
        if(vparent <= vchild),
            basis(node(d,b)) = 0;
            value(node(d,b)) = vparent;
        else
            basis(node(d,b)) = 1;
            value(node(d,b)) = vchild;
        end
    end
end
if strcmp(WLVERBOSE,'Yes')
    fprintf('best basis %g \n',value(1))
end


%
% Copyright (c) 1993. David L. Donoho
%


%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%
function index = node(d,b)
% node -- Tree indexing function
%  Usage
%    index = node(d,b)
%  Inputs
%    d        depth from root of tree
%    b        index among the 2^d possibilities
%             in a left-right scan at that depth
%  Outputs
%    index    linear index of node in tree structure
%
	index =  2^d + b;

%
% Copyright (c) 1993. David L. Donoho
%     
    
    
    
 
 
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

xi = (1 + (.5:(m-.5))./m)./2;
if strcmp(Name,'Trivial'),
    bp = sqrt(xi);
elseif strcmp(Name,'Sine'),
    bp = sin( pi/2 .* xi );
end
bm = sqrt(1 - bp .^2);

%
% Copyright (c) 1993. David L. Donoho
%


%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%



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
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%



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
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%



function x = ipt_cp(basis,coef,D,bellname)
% IPT_CP -- Synthesize signal from cosine packet coefficients
%  Usage
%    x = IPT_CP(btree,coef,D,bell)
%  Inputs
%    btree    basis tree selecting Cosine Packet basis
%    coef     coefficients in that basis
%    D        maximum splitting depth
%    bell     name of orthonormal bell
%  Outputs
%    x        1-d signal whose cosine packet coeff's in
%             basis btree come from cp
%
%  Description
%    Perform the inverse operation of FPT_CP.
%
%  See Also
%      CPAnalysis, CPTour, MakeONBell
%
[n,J] = dyadlength(coef);

% Create Bell
if nargin < 3,
    bellname = 'Sine';
end
m = n / 2^D /2;
[bp,bm] = MakeONBell(bellname,m);

% Setup arrays as rows
x  = zeros(1,n);
cf = ShapeAsRow(coef);

% initialize tree traversal stack
stack = zeros(2,2^D+1);
k = 1;
stack(:,k) = [0 0 ]';

while(k > 0),

    % pop stack
    d = stack(1,k);
    b = stack(2,k);
    k = k-1;

    if(basis(node(d,b)) ~= 0) ,  % nonterminal node
        k = k+1; stack(:,k) = [(d+1) (2*b)  ]';
        k = k+1; stack(:,k) = [(d+1) (2*b+1)]';
    else
        c = cf(packet(d,b,n));
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
end
x = ShapeLike(x,coef);

%
% Copyright (c) 1993. David L. Donoho
%


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
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%

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


function vec = ShapeLike(sig,proto)
% ShapeLike -- Make 1-d signal with given shape
%  Usage
%    vec = ShapeLike(sig,proto)
%  Inputs
%    sig      a row or column vector
%    proto    a prototype shape (row or column vector)
%  Outputs
%    vec      a vector with contents taken from sig
%             and same shape as proto
%
%  See Also
%    ShapeAsRow
%
	sp = size(proto);
	ss = size(sig);
	if( sp(1)>1 & sp(2)>1 )
	   disp('Weird proto argument to ShapeLike')
	elseif ss(1)>1 & ss(2) > 1,
	   disp('Weird sig argument to ShapeLike')
	else
	   if(sp(1) > 1),
		  if ss(1) > 1,
			 vec = sig;
		  else
			 vec = sig(:);
		  end
	   else
		  if ss(2) > 1,
			 vec = sig;
		  else
			 vec = sig(:)';
		  end
	   end
	end
    
    
 
 
%
%  Part of Wavelab Version 850
%  Built Tue Jan  3 13:20:43 EST 2006
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
