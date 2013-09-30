function y = perform_wavelet_transform(x, Jmin, dir, options)

% perform_wavelet_transform - wrapper to wavelab Wavelet transform (1D/2D and orthogonal/biorthogonal).
%
%   y = perform_wavelet_transform(x, Jmin, dir, options);
%
%   'x' is either a 1D or a 2D array.
%   'Jmin' is the minimum scale (i.e. the coarse channel is of size 2^Jmin
%       in 1D).
%   'dir' is +1 for fwd transform and -1 for bwd.
%   'options.wavelet_vm' is the number of Vanishing moment (both for primal and dual).
%   'options.wavelet_type' can be
%       'daubechies', 'symmlet', 'battle', 'biorthogonal'.
%
%   Typical use : 
%       M = <load your image here>;
%       Jmin = 4;
%       options.wavelet_type = 'biorthogonal_swapped';
%       options.wavelet_vm = 4;
%       MW = perform_wavelet_transform(M, Jmin, +1, options);
%       Mt = <perform some modification on MW>
%       M = perform_wavelet_transform(Mt, Jmin, -1, options);
%
%   'y' is an array of the same size as 'x'. This means that for the 2D
%   we are stuck to the wavelab coding style, i.e. the result
%   of each transform is an array organized using Mallat's ordering
%   (whereas Matlab official toolbox use a 1D ordering for the 2D transform).
%
%   Here the transform automaticaly select symmetric boundary condition
%   if you use a symmetric filter. If your filter is not symmetric
%   (e.g. Dauechies filters) then as the output must have same length
%   as the input, the boundary condition are automatically set to periodic.
%
%   You do not need Wavelab to use this function (the Wavelab .m file are
%   included in this script). However, for faster execution time, you
%   should install the mex file within the Wavelab distribution.
%       http://www-stat.stanford.edu/~wavelab/
%
%   Copyright (c) 2005 Gabriel Peyre

if nargin<3
    dir = 1;
end
if nargin<2
    Jmin = 3;
end

options.null = 0;

if isfield(options, 'wavelet_type')
    wavelet_type = options.wavelet_type;
else
    wavelet_type = 'daubechies';
end

if isfield(options, 'wavelet_vm')
    VM = options.wavelet_vm;
else
    VM = 4;
end


if isfield(options, 'ti')
    % for translation-invariant transform
    ti = options.ti;
else
    ti = 0;
end

ndim = length(size(x));
if ndim==2 && ( size(x,2)==1 || size(x,1)==1 )
    ndim=1;
end

% for color images
if ndims(x)>2
    y = x;
    for i=1:size(x,3)
        y(:,:,i) = perform_wavelet_transform(x(:,:,i), Jmin, dir, options);
    end
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           WAVELAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate filters
switch lower(wavelet_type)
    case 'daubechies'
        qmf = MakeONFilter('Daubechies',VM*2);  % in Wavelab, 2nd argument is VM*2 for Daubechies... no comment ...
    case 'haar'
        qmf = MakeONFilter('Haar');  % in Wavelab, 2nd argument is VM*2 for Daubechies... no comment ...
    case 'symmlet'
        qmf = MakeONFilter('Symmlet',VM);
    case 'battle'
        qmf = MakeONFilter('Battle',VM-1);
        dqmf = qmf; % we need dual filter
    case 'biorthogonal'
        [qmf,dqmf] = MakeBSFilter( 'CDF', [VM,VM] );
    case 'biorthogonal_swapped'
        [dqmf,qmf] = MakeBSFilter( 'CDF', [VM,VM] );
    otherwise
        error('Unknown transform.');
end

% translation invariant transform
if ti
    if ndim==2 && size(x,2)<50
        ndim = 1;
    end
    if ndim==1
        if dir==1
            y = TI2Stat( FWT_TI(x,Jmin,qmf) );
        else
            y = IWT_TI( Stat2TI(x),qmf);
        end
    elseif ndim==2
        if dir==1
            y = FWT2_TI(x,Jmin,qmf);
        else
            y = IWT2_TI(x,Jmin,qmf);
        end
    end
    return;
end

% perform transform
if ~exist('dqmf')
    %%% ORTHOGONAL %%%
    if ndim==1
        if dir==1
            y = FWT_PO(x,Jmin,qmf);
        else
            y = IWT_PO(x,Jmin,qmf);
        end
    elseif ndim==2
        if dir==1
            y = FWT2_PO(x,Jmin,qmf);
        else
            y = IWT2_PO(x,Jmin,qmf);
        end
    end
else
    %%% BI-ORTHOGONAL %%%
    if ndim==1
        if dir==1
            y = FWT_SBS(x,Jmin,qmf,dqmf);
        else
            y = IWT_SBS(x,Jmin,qmf,dqmf);
        end
    elseif ndim==2
        if dir==1
            y = FWT2_SBS(x,Jmin,qmf,dqmf);
        else
            y = IWT2_SBS(x,Jmin,qmf,dqmf);
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WAVELAB DISTRIBUTION -- http://www-stat.stanford.edu/~wavelab/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WAVELAB DISTRIBUTION -- http://www-stat.stanford.edu/~wavelab/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function f = MakeONFilter(Type,Par)
% MakeONFilter -- Generate Orthonormal QMF Filter for Wavelet Transform
%  Usage
%    qmf = MakeONFilter(Type,Par)
%  Inputs
%    Type   string, 'Haar', 'Beylkin', 'Coiflet', 'Daubechies',
%           'Symmlet', 'Vaidyanathan','Battle'
%    Par    integer, it is a parameter related to the support and vanishing
%           moments of the wavelets, explained below for each wavelet.
%
% Outputs
%    qmf    quadrature mirror filter
%
%  Description
%    The Haar filter (which could be considered a Daubechies-2) was the
%    first wavelet, though not called as such, and is discontinuous.
%
%    The Beylkin filter places roots for the frequency response function
%    close to the Nyquist frequency on the real axis.
%
%    The Coiflet filters are designed to give both the mother and father
%    wavelets 2*Par vanishing moments; here Par may be one of 1,2,3,4 or 5.
%
%    The Daubechies filters are minimal phase filters that generate wavelets
%    which have a minimal support for a given number of vanishing moments.
%    They are indexed by their length, Par, which may be one of
%    4,6,8,10,12,14,16,18 or 20. The number of vanishing moments is par/2.
%
%    Symmlets are also wavelets within a minimum size support for a given
%    number of vanishing moments, but they are as symmetrical as possible,
%    as opposed to the Daubechies filters which are highly asymmetrical.
%    They are indexed by Par, which specifies the number of vanishing
%    moments and is equal to half the size of the support. It ranges
%    from 4 to 10.
%
%    The Vaidyanathan filter gives an exact reconstruction, but does not
%    satisfy any moment condition.  The filter has been optimized for
%    speech coding.
%
%    The Battle-Lemarie filter generate spline orthogonal wavelet basis.
%    The parameter Par gives the degree of the spline. The number of
%    vanishing moments is Par+1.
%
%  See Also
%    FWT_PO, IWT_PO, FWT2_PO, IWT2_PO, WPAnalysis
%
%  References
%    The books by Daubechies and Wickerhauser.
%

if strcmp(Type,'Haar'),
    f = [1 1] ./ sqrt(2);
end

if strcmp(Type,'Beylkin'),
    f = [	.099305765374	.424215360813	.699825214057	...
        .449718251149	-.110927598348	-.264497231446	...
        .026900308804	.155538731877	-.017520746267	...
        -.088543630623	.019679866044	.042916387274	...
        -.017460408696	-.014365807969	.010040411845	...
        .001484234782	-.002736031626	.000640485329	];
end

if strcmp(Type,'Coiflet'),
    if Par==1,
        f = [	.038580777748	-.126969125396	-.077161555496	...
            .607491641386	.745687558934	.226584265197	];
    end
    if Par==2,
        f = [	.016387336463	-.041464936782	-.067372554722	...
            .386110066823	.812723635450	.417005184424	...
            -.076488599078	-.059434418646	.023680171947	...
            .005611434819	-.001823208871	-.000720549445	];
    end
    if Par==3,
        f = [	-.003793512864	.007782596426	.023452696142	...
            -.065771911281	-.061123390003	.405176902410	...
            .793777222626	.428483476378	-.071799821619	...
            -.082301927106	.034555027573	.015880544864	...
            -.009007976137	-.002574517688	.001117518771	...
            .000466216960	-.000070983303	-.000034599773	];
    end
    if Par==4,
        f = [	.000892313668	-.001629492013	-.007346166328	...
            .016068943964	.026682300156	-.081266699680	...
            -.056077313316	.415308407030	.782238930920	...
            .434386056491	-.066627474263	-.096220442034	...
            .039334427123	.025082261845	-.015211731527	...
            -.005658286686	.003751436157	.001266561929	...
            -.000589020757	-.000259974552	.000062339034	...
            .000031229876	-.000003259680	-.000001784985	];
    end
    if Par==5,
        f = [	-.000212080863	.000358589677	.002178236305	...
            -.004159358782	-.010131117538	.023408156762	...
            .028168029062	-.091920010549	-.052043163216	...
            .421566206729	.774289603740	.437991626228	...
            -.062035963906	-.105574208706	.041289208741	...
            .032683574283	-.019761779012	-.009164231153	...
            .006764185419	.002433373209	-.001662863769	...
            -.000638131296	.000302259520	.000140541149	...
            -.000041340484	-.000021315014	.000003734597	...
            .000002063806	-.000000167408	-.000000095158	];
    end
end

if strcmp(Type,'Daubechies'),
    if Par==4,
        f = [	.482962913145	.836516303738	...
            .224143868042	-.129409522551	];
    end
    if Par==6,
        f = [	.332670552950	.806891509311	...
            .459877502118	-.135011020010	...
            -.085441273882	.035226291882	];
    end
    if Par==8,
        f = [ 	.230377813309	.714846570553	...
            .630880767930	-.027983769417	...
            -.187034811719	.030841381836	...
            .032883011667	-.010597401785	];
    end
    if Par==10,
        f = [	.160102397974	.603829269797	.724308528438	...
            .138428145901	-.242294887066	-.032244869585	...
            .077571493840	-.006241490213	-.012580751999	...
            .003335725285									];
    end
    if Par==12,
        f = [	.111540743350	.494623890398	.751133908021	...
            .315250351709	-.226264693965	-.129766867567	...
            .097501605587	.027522865530	-.031582039317	...
            .000553842201	.004777257511	-.001077301085	];
    end
    if Par==14,
        f = [	.077852054085	.396539319482	.729132090846	...
            .469782287405	-.143906003929	-.224036184994	...
            .071309219267	.080612609151	-.038029936935	...
            -.016574541631	.012550998556	.000429577973	...
            -.001801640704	.000353713800					];
    end
    if Par==16,
        f = [	.054415842243	.312871590914	.675630736297	...
            .585354683654	-.015829105256	-.284015542962	...
            .000472484574	.128747426620	-.017369301002	...
            -.044088253931	.013981027917	.008746094047	...
            -.004870352993	-.000391740373	.000675449406	...
            -.000117476784									];
    end
    if Par==18,
        f = [	.038077947364	.243834674613	.604823123690	...
            .657288078051	.133197385825	-.293273783279	...
            -.096840783223	.148540749338	.030725681479	...
            -.067632829061	.000250947115	.022361662124	...
            -.004723204758	-.004281503682	.001847646883	...
            .000230385764	-.000251963189	.000039347320	];
    end
    if Par==20,
        f = [	.026670057901	.188176800078	.527201188932	...
            .688459039454	.281172343661	-.249846424327	...
            -.195946274377	.127369340336	.093057364604	...
            -.071394147166	-.029457536822	.033212674059	...
            .003606553567	-.010733175483	.001395351747	...
            .001992405295	-.000685856695	-.000116466855	...
            .000093588670	-.000013264203					];
    end
end

if strcmp(Type,'Symmlet'),
    if Par==4,
        f = [	-.107148901418	-.041910965125	.703739068656	...
            1.136658243408	.421234534204	-.140317624179	...
            -.017824701442	.045570345896					];
    end
    if Par==5,
        f = [	.038654795955	.041746864422	-.055344186117	...
            .281990696854	1.023052966894	.896581648380	...
            .023478923136	-.247951362613	-.029842499869	...
            .027632152958									];
    end
    if Par==6,
        f = [	.021784700327	.004936612372	-.166863215412	...
            -.068323121587	.694457972958	1.113892783926	...
            .477904371333	-.102724969862	-.029783751299	...
            .063250562660	.002499922093	-.011031867509	];
    end
    if Par==7,
        f = [	.003792658534	-.001481225915	-.017870431651	...
            .043155452582	.096014767936	-.070078291222	...
            .024665659489	.758162601964	1.085782709814	...
            .408183939725	-.198056706807	-.152463871896	...
            .005671342686	.014521394762					];
    end
    if Par==8,
        f = [	.002672793393	-.000428394300	-.021145686528	...
            .005386388754	.069490465911	-.038493521263	...
            -.073462508761	.515398670374	1.099106630537	...
            .680745347190	-.086653615406	-.202648655286	...
            .010758611751	.044823623042	-.000766690896	...
            -.004783458512									];
    end
    if Par==9,
        f = [	.001512487309	-.000669141509	-.014515578553	...
            .012528896242	.087791251554	-.025786445930	...
            -.270893783503	.049882830959	.873048407349	...
            1.015259790832	.337658923602	-.077172161097	...
            .000825140929	.042744433602	-.016303351226	...
            -.018769396836	.000876502539	.001981193736	];
    end
    if Par==10,
        f = [	.001089170447	.000135245020	-.012220642630	...
            -.002072363923	.064950924579	.016418869426	...
            -.225558972234	-.100240215031	.667071338154	...
            1.088251530500	.542813011213	-.050256540092	...
            -.045240772218	.070703567550	.008152816799	...
            -.028786231926	-.001137535314	.006495728375	...
            .000080661204	-.000649589896					];
    end
end

if strcmp(Type,'Vaidyanathan'),
    f = [	-.000062906118	.000343631905	-.000453956620	...
        -.000944897136	.002843834547	.000708137504	...
        -.008839103409	.003153847056	.019687215010	...
        -.014853448005	-.035470398607	.038742619293	...
        .055892523691	-.077709750902	-.083928884366	...
        .131971661417	.135084227129	-.194450471766	...
        -.263494802488	.201612161775	.635601059872	...
        .572797793211	.250184129505	.045799334111		];
end

if strcmp(Type,'Battle'),
    if Par == 1,
        g = [0.578163    0.280931   -0.0488618   -0.0367309 ...
            0.012003    0.00706442 -0.00274588 -0.00155701 ...
            0.000652922 0.000361781 -0.000158601 -0.0000867523
            ];
    end

    if Par == 3,

        g = [0.541736    0.30683    -0.035498    -0.0778079 ...
            0.0226846   0.0297468     -0.0121455 -0.0127154 ...
            0.00614143 0.00579932    -0.00307863 -0.00274529 ...
            0.00154624 0.00133086 -0.000780468 -0.00065562 ...
            0.000395946 0.000326749 -0.000201818 -0.000164264 ...
            0.000103307
            ];
    end

    if Par == 5,
        g = [0.528374    0.312869    -0.0261771   -0.0914068 ...
            0.0208414    0.0433544 -0.0148537 -0.0229951  ...
            0.00990635 0.0128754    -0.00639886 -0.00746848 ...
            0.00407882 0.00444002 -0.00258816    -0.00268646 ...
            0.00164132 0.00164659 -0.00104207 -0.00101912 ...
            0.000662836 0.000635563 -0.000422485 -0.000398759 ...
            0.000269842 0.000251419 -0.000172685 -0.000159168 ...
            0.000110709 0.000101113
            ];
    end
    l = length(g);
    f = zeros(1,2*l-1);
    f(l:2*l-1) = g;
    f(1:l-1) = reverse(g(2:l));
end

f = f ./ norm(f);

%
% Copyright (c) 1993-5. Jonathan Buckheit and David Donoho
%


%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%

function wcoef = FWT_PO(x,L,qmf)
% FWT_PO -- Forward Wavelet Transform (periodized, orthogonal)
%  Usage
%    wc = FWT_PO(x,L,qmf)
%  Inputs
%    x    1-d signal; length(x) = 2^J
%    L    Coarsest Level of V_0;  L << J
%    qmf  quadrature mirror filter (orthonormal)
%  Outputs
%    wc    1-d wavelet transform of x.
%
%  Description
%    1. qmf filter may be obtained from MakeONFilter
%    2. usually, length(qmf) < 2^(L+1)
%    3. To reconstruct use IWT_PO
%
%  See Also
%    IWT_PO, MakeONFilter
%
[n,J] = dyadlength(x) ;
wcoef = zeros(1,n) ;
beta = ShapeAsRow(x);  %take samples at finest scale as beta-coeffts
for j=J-1:-1:L
    alfa = DownDyadHi(beta,qmf);
    wcoef(dyad(j)) = alfa;
    beta = DownDyadLo(beta,qmf) ;
end
wcoef(1:(2^L)) = beta;
wcoef = ShapeLike(wcoef,x);

%
% Copyright (c) 1993. Iain M. Johnstone
%


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

function d = DownDyadHi(x,qmf)
% DownDyadHi -- Hi-Pass Downsampling operator (periodized)
%  Usage
%    d = DownDyadHi(x,f)
%  Inputs
%    x    1-d signal at fine scale
%    f    filter
%  Outputs
%    y    1-d signal at coarse scale
%
%  See Also
%    DownDyadLo, UpDyadHi, UpDyadLo, FWT_PO, iconv
%
d = iconv( MirrorFilt(qmf),lshift(x));
n = length(d);
d = d(1:2:(n-1));

%
% Copyright (c) 1993. Iain M. Johnstone
%


%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%


function y = MirrorFilt(x)
% MirrorFilt -- Apply (-1)^t modulation
%  Usage
%    h = MirrorFilt(l)
%  Inputs
%    l   1-d signal
%  Outputs
%    h   1-d signal with DC frequency content shifted
%        to Nyquist frequency
%
%  Description
%    h(t) = (-1)^(t-1)  * x(t),  1 <= t <= length(x)
%
%  See Also
%    DyadDownHi
%

y = -( (-1).^(1:length(x)) ).*x;

%
% Copyright (c) 1993. Iain M. Johnstone
%


%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%


function y = lshift(x)
% lshift -- Circular left shift of 1-d signal
%  Usage
%    l = lshift(x)
%  Inputs
%    x   1-d signal
%  Outputs
%    l   1-d signal
%        l(i) = x(i+1) except l(n) = x(1)
%

y = [ x( 2:length(x) ) x(1) ];

%
% Copyright (c) 1993. Iain M. Johnstone
%


%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%


function y = iconv(f,x)
% iconv -- Convolution Tool for Two-Scale Transform
%  Usage
%    y = iconv(f,x)
%  Inputs
%    f   filter
%    x   1-d signal
%  Outputs
%    y   filtered result
%
%  Description
%    Filtering by periodic convolution of x with f
%
%  See Also
%    aconv, UpDyadHi, UpDyadLo, DownDyadHi, DownDyadLo
%
n = length(x);
p = length(f);
if p <= n,
    xpadded = [x((n+1-p):n) x];
else
    z = zeros(1,p);
    for i=1:p,
        imod = 1 + rem(p*n -p + i-1,n);
        z(i) = x(imod);
    end
    xpadded = [z x];
end
ypadded = filter(f,1,xpadded);
y = ypadded((p+1):(n+p));

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

function i = dyad(j)
% dyad -- Index entire j-th dyad of 1-d wavelet xform
%  Usage
%    ix = dyad(j);
%  Inputs
%    j     integer
%  Outputs
%    ix    list of all indices of wavelet coeffts at j-th level
%
i = (2^(j)+1):(2^(j+1)) ;

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

function d = DownDyadLo(x,qmf)
% DownDyadLo -- Lo-Pass Downsampling operator (periodized)
%  Usage
%    d = DownDyadLo(x,f)
%  Inputs
%    x    1-d signal at fine scale
%    f    filter
%  Outputs
%    y    1-d signal at coarse scale
%
%  See Also
%    DownDyadHi, UpDyadHi, UpDyadLo, FWT_PO, aconv
%
d = aconv(qmf,x);
n = length(d);
d = d(1:2:(n-1));

%
% Copyright (c) 1993. Iain M. Johnstone
%


%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%

function y = aconv(f,x)
% aconv -- Convolution Tool for Two-Scale Transform
%  Usage
%    y = aconv(f,x)
%  Inputs
%    f    filter
%    x    1-d signal
%  Outputs
%    y    filtered result
%
%  Description
%    Filtering by periodic convolution of x with the
%    time-reverse of f.
%
%  See Also
%    iconv, UpDyadHi, UpDyadLo, DownDyadHi, DownDyadLo
%

n = length(x);
p = length(f);
if p < n,
    xpadded = [x x(1:p)];
else
    z = zeros(1,p);
    for i=1:p,
        imod = 1 + rem(i-1,n);
        z(i) = x(imod);
    end
    xpadded = [x z];
end
fflip = reverse(f);
ypadded = filter(fflip,1,xpadded);
y = ypadded(p:(n+p-1));

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
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%

function x = IWT_PO(wc,L,qmf)
% IWT_PO -- Inverse Wavelet Transform (periodized, orthogonal)
%  Usage
%    x = IWT_PO(wc,L,qmf)
%  Inputs
%    wc     1-d wavelet transform: length(wc) = 2^J.
%    L      Coarsest scale (2^(-L) = scale of V_0); L << J;
%    qmf    quadrature mirror filter
%  Outputs
%    x      1-d signal reconstructed from wc
%
%  Description
%    Suppose wc = FWT_PO(x,L,qmf) where qmf is an orthonormal quad. mirror
%    filter, e.g. one made by MakeONFilter. Then x can be reconstructed by
%      x = IWT_PO(wc,L,qmf)
%
%  See Also
%    FWT_PO, MakeONFilter
%
wcoef = ShapeAsRow(wc);
x = wcoef(1:2^L);
[n,J] = dyadlength(wcoef);
for j=L:J-1
    x = UpDyadLo(x,qmf) + UpDyadHi(wcoef(dyad(j)),qmf)  ;
end
x = ShapeLike(x,wc);

%
% Copyright (c) 1993. Iain M. Johnstone
%


%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%

function y = UpDyadLo(x,qmf)
% UpDyadLo -- Lo-Pass Upsampling operator; periodized
%  Usage
%    u = UpDyadLo(d,f)
%  Inputs
%    d    1-d signal at coarser scale
%    f    filter
%  Outputs
%    u    1-d signal at finer scale
%
%  See Also
%    DownDyadLo, DownDyadHi, UpDyadHi, IWT_PO, iconv
%
y =  iconv(qmf, UpSample(x) );

%
% Copyright (c) 1993. Iain M. Johnstone
%


%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%

function y = UpSample(x,s)
% UpSample -- Upsampling operator
%  Usage
%    u = UpSample(d[,s])
%  Inputs
%    d   1-d signal, of length n
%    s   upsampling scale, default = 2
%  Outputs
%    u   1-d signal, of length s*n with zeros
%        interpolating alternate samples
%        u(s*i-1) = d(i), i=1,...,n
%

if nargin == 1, s = 2; end
n = length(x)*s;
y = zeros(1,n);
y(1:s:(n-s+1) )=x;


%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%

function y = UpDyadHi(x,qmf)
% UpDyadHi -- Hi-Pass Upsampling operator; periodized
%  Usage
%    u = UpDyadHi(d,f)
%  Inputs
%    d    1-d signal at coarser scale
%    f    filter
%  Outputs
%    u    1-d signal at finer scale
%
%  See Also
%    DownDyadLo, DownDyadHi, UpDyadLo, IWT_PO, aconv
%

y = aconv( MirrorFilt(qmf), rshift( UpSample(x) ) );

%
% Copyright (c) 1993. Iain M. Johnstone
%


%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%

function y = rshift(x)
% rshift -- Circular right shift of 1-d signal
%  Usage
%    r = rshift(x)
%  Inputs
%    x   1-d signal
%  Outputs
%    r   1-d signal
%        r(i) = x(i-1) except r(1) = x(n)
%

n = length(x);
y = [ x(n) x( 1: (n-1) )];

%
% Copyright (c) 1993. Iain M. Johnstone
%


%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%




function [qmf,dqmf] = MakeBSFilter(Type,Par)
% MakeBSFilter -- Generate Biorthonormal QMF Filter Pairs
%  Usage
%    [qmf,dqmf] = MakeBSFilter(Type,Par)
%  Inputs
%    Type   string, one of:
%             'Triangle'
%             'Interpolating' 'Deslauriers' (the two are same)
%             'Average-Interpolating'
%             'CDF' (spline biorthogonal filters in Daubechies's book)
%             'Villasenor' (Villasenor's 5 best biorthogonal filters)
%    Par    integer list, e.g. if Type ='Deslauriers', Par=3 specifies
%           Deslauriers-Dubuc filter, polynomial degree 3
%  Outputs
%    qmf    quadrature mirror filter  (odd length, symmetric)
%    dqmf   dual quadrature mirror filter  (odd length, symmetric)
%
%  See Also
%    FWT_PBS, IWT_PBS, FWT2_PBS, IWT2_PBS
%
%  References
%    I. Daubechies, "Ten Lectures on Wavelets."
%
%    G. Deslauriers and S. Dubuc, "Symmetric Iterative Interpolating Processes."
%
%    D. Donoho, "Smooth Wavelet Decompositions with Blocky Coefficient Kernels."
%
%    J. Villasenor, B. Belzer and J. Liao, "Wavelet Filter Evaluation for
%    Image Compression."
%

if nargin < 2,
    Par = 0;
end

sqr2 = sqrt(2);

if strcmp(Type,'Triangle'),
    qmf = [0 1 0];
    dqmf = [.5 1 .5];

elseif strcmp(Type,'Interpolating') | strcmp(Type,'Deslauriers'),
    qmf  = [0 1 0];
    dqmf = MakeDDFilter(Par)';
    dqmf =  dqmf(1:(length(dqmf)-1));

elseif strcmp(Type,'Average-Interpolating'),
    qmf  = [0 .5 .5] ;
    dqmf = [0 ; MakeAIFilter(Par)]';

elseif strcmp(Type,'CDF'),
    if Par(1)==1,
        dqmf = [0 .5 .5] .* sqr2;
        if Par(2) == 1,
            qmf = [0 .5 .5] .* sqr2;
        elseif Par(2) == 3,
            qmf = [0 -1 1 8 8 1 -1] .* sqr2 / 16;
        elseif Par(2) == 5,
            qmf = [0 3 -3 -22 22 128 128 22 -22 -3 3].*sqr2/256;
        end
    elseif Par(1)==2,
        dqmf = [.25 .5 .25] .* sqr2;
        if Par(2)==2,
            qmf = [-.125 .25 .75 .25 -.125] .* sqr2;
        elseif Par(2)==4,
            qmf = [3 -6 -16 38 90 38 -16 -6 3] .* (sqr2/128);
        elseif Par(2)==6,
            qmf = [-5 10 34 -78 -123 324 700 324 -123 -78 34 10 -5 ] .* (sqr2/1024);
        elseif Par(2)==8,
            qmf = [35 -70 -300 670 1228 -3126 -3796 10718 22050 ...
                10718 -3796 -3126 1228 670 -300 -70 35 ] .* (sqr2/32768);
        end
    elseif Par(1)==3,
        dqmf = [0 .125 .375 .375 .125] .* sqr2;
        if Par(2) == 1,
            qmf = [0 -.25 .75 .75 -.25] .* sqr2;
        elseif Par(2) == 3,
            qmf = [0 3 -9 -7 45 45 -7 -9 3] .* sqr2/64;
        elseif Par(2) == 5,
            qmf = [0 -5 15 19 -97 -26 350 350 -26 -97 19 15 -5] .* sqr2/512;
        elseif Par(2) == 7,
            qmf = [0 35 -105 -195 865 363 -3489 -307 11025 11025 -307 -3489 363 865 -195 -105 35] .* sqr2/16384;
        elseif Par(2) == 9,
            qmf = [0 -63 189 469 -1911 -1308 9188 1140 -29676 190 87318 87318 190 -29676 ...
                1140 9188 -1308 -1911 469 189 -63] .* sqr2/131072;
        end
    elseif Par(1)==4,
        dqmf = [.026748757411 -.016864118443 -.078223266529 .266864118443 .602949018236 ...
            .266864118443 -.078223266529 -.016864118443 .026748757411] .*sqr2;
        if Par(2) == 4,
            qmf = [0 -.045635881557 -.028771763114 .295635881557 .557543526229 ...
                .295635881557 -.028771763114 -.045635881557 0] .*sqr2;
        end
    end

elseif strcmp(Type,'Villasenor'),
    if Par == 1,
        % The "7-9 filters"
        qmf = [.037828455506995 -.023849465019380 -.11062440441842 .37740285561265];
        qmf = [qmf .85269867900940 reverse(qmf)];
        dqmf = [-.064538882628938 -.040689417609558 .41809227322221];
        dqmf = [dqmf .78848561640566 reverse(dqmf)];
    elseif Par == 2,
        qmf  = [-.008473 .003759 .047282 -.033475 -.068867 .383269 .767245 .383269 -.068867...
            -.033475 .047282 .003759 -.008473];
        dqmf = [0.014182  0.006292 -0.108737 -0.069163 0.448109 .832848 .448109 -.069163 -.108737 .006292 .014182];
    elseif Par == 3,
        qmf  = [0 -.129078 .047699 .788486 .788486 .047699 -.129078];
        dqmf = [0 .018914 .006989 -.067237 .133389 .615051 .615051 .133389 -.067237 .006989 .018914];
    elseif Par == 4,
        qmf  = [-1 2 6 2 -1] / (4*sqr2);
        dqmf = [1 2 1] / (2*sqr2);
    elseif Par == 5,
        qmf  = [0 1 1]/sqr2;
        dqmf = [0 -1 1 8 8 1 -1]/(8*sqr2);
    end
end

%
% Copyright (c) 1995. Jonathan Buckheit, Shaobing Chen and David Donoho
%
% Modified by Maureen Clerc and Jerome Kalifa, 1997
% clerc@cmapx.polytechnique.fr, kalifa@cmapx.polytechnique.fr
%


%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%

function wcoef = FWT_SBS(x,L,qmf,dqmf)
% FWT_SBS -- Forward Wavelet Transform (symmetric extension, biorthogonal, symmetric)
%  Usage
%    wc = FWT_SBS(x,L,qmf,dqmf)
%  Inputs
%    x    1-d signal; arbitrary length
%    L    Coarsest Level of V_0;  L << J
%    qmf    quadrature mirror filter (symmetric)
%    dqmf   quadrature mirror filter (symmetric, dual of qmf)
%  Outputs
%    wc    1-d wavelet transform of x.
%
%  Description
%    1. qmf filter may be obtained from MakePBSFilter
%    2. usually, length(qmf) < 2^(L+1)
%    3. To reconstruct use IWT_SBS
%
%  See Also
%    IWT_SBS, MakePBSFilter
%
%  References
%   Based on the algorithm developed by Christopher Brislawn.
%   See "Classification of Symmetric Wavelet Transforms"
%

[n,J] = dyadlength(x);

wcoef = zeros(1,n);
beta = ShapeAsRow(x);  % take samples at finest scale as beta-coeffts

dp = dyadpartition(n);

for j=J-1:-1:L,
    [beta, alfa] = DownDyad_SBS(beta,qmf,dqmf);
    dyadj = (dp(j+1)+1):dp(j+2);
    wcoef(dyadj) = alfa;
end
wcoef(1:length(beta)) = beta;
wcoef = ShapeLike(wcoef,x);

%
% Copyright (c) 1996. Thomas P.Y. Yu
%



%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%

function dp = dyadpartition(n)
% dyadpartition -- determine dyadic partition in wavelet transform of
%                  nondyadic signals

J = ceil(log2(n));

m = n;
for j=J-1:-1:0;
    if rem(m,2)==0,
        dps(j+1) = m/2;
        m = m/2;
    else
        dps(j+1) = (m-1)/2;
        m = (m+1)/2;
    end
end

dp = cumsum([1 dps]);

%
% Copyright (c) 1996. Thomas P.Y. Yu
%


%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%

function [beta, alpha] = DownDyad_SBS(x,qmf,dqmf)
% DownDyad_SBS -- Symmetric Downsampling operator
%  Usage
%    [beta,alpha] = DownDyad_SBS(x,qmf,dqmf)
%  Inputs
%    x     1-d signal at fine scale
%    qmf   quadrature mirror filter
%    dqmf  dual quadrature mirror filter
%  Outputs
%    beta  coarse coefficients
%    alpha fine coefficients
%  See Also
%    FWT_SBS
%

% oddf = (rem(length(qmf),2)==1);
oddf = ~(qmf(1)==0 & qmf(length(qmf))~=0);
oddx = (rem(length(x),2)==1);

% symmetric extension of x
if oddf,
    ex = extend(x,1,1);
else
    ex = extend(x,2,2);
end

% convolution
ebeta = DownDyadLo_PBS(ex,qmf);
ealpha = DownDyadHi_PBS(ex,dqmf);

% project
if oddx,
    beta = ebeta(1:(length(x)+1)/2);
    alpha = ealpha(1:(length(x)-1)/2);
else
    beta = ebeta(1:length(x)/2);
    alpha = ealpha(1:length(x)/2);
end

%
% Copyright (c) 1996. Thomas P.Y. Yu
%


%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%

function y = extend(x, par1, par2)
% extend -- perform various kinds of symmetric extension
%
if par1==1 & par2==1,
    y = [x x((length(x)-1):-1:2)];
elseif par1==1 & par2==2,
    y = [x x((length(x)-1):-1:1)];
elseif par1==2 & par2==1,
    y = [x x(length(x):-1:2)];
elseif par1==2 & par2==2,
    y = [x reverse(x)];
end

%
% Copyright (c) 1996. Thomas P.Y. Yu
%


%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%

function d = DownDyadLo_PBS(x,qmf)
% DownDyadLo_PBS -- Lo-Pass Downsampling operator (periodized,symmetric)
%  Usage
%    d = DownDyadLo_PBS(x,sf)
%  Inputs
%    x    1-d signal at fine scale
%    sf   symmetric filter
%  Outputs
%    y    1-d signal at coarse scale
%
%  See Also
%    DownDyadHi_PBS, UpDyadHi_PBS, UpDyadLo_PBS, FWT_PBSi, symm_aconv
%
d = symm_aconv(qmf,x);
n = length(d);
d = d(1:2:(n-1));

%
% Copyright (c) 1995. David L. Donoho
%




%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%

function y = symm_aconv(sf,x)
% symm_aconv -- Symmetric Convolution Tool for Two-Scale Transform
%  Usage
%    y = symm_aconv(sf,x)
%  Inputs
%    sf   symmetric filter
%    x    1-d signal
%  Outputs
%    y    filtered result
%
%  Description
%    Filtering by periodic convolution of x with the
%    time-reverse of sf.
%
%  See Also
%    symm_iconv, UpDyadHi_PBS, UpDyadLo_PBS, DownDyadHi_PBS, DownDyadLo_PBS
%

n = length(x);
p = length(sf);
if p < n,
    xpadded = [x x(1:p)];
else
    z = zeros(1,p);
    for i=1:p,
        imod = 1 + rem(i-1,n);
        z(i) = x(imod);
    end
    xpadded = [x z];
end

fflip = reverse(sf);
ypadded = filter(fflip,1,xpadded);

if p < n,
    y = [ypadded((n+1):(n+p)) ypadded((p+1):(n))];
else
    for i=1:n,
        imod = 1 + rem(p+i-1,n);
        y(imod) = ypadded(p+i);
    end
end

shift = (p-1)/ 2 ;
shift = 1 + rem(shift-1, n);
y = [y((1+shift):n) y(1:(shift))] ;

%
% Copyright (c) 1995. Shaobing Chen and David L. Donoho
%




%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%

function d = DownDyadHi_PBS(x,qmf)
% DownDyadHi_PBS -- Hi-Pass Downsampling operator (periodized,symmetric)
%  Usage
%    d = DownDyadHi_PBS(x,sqmf)
%  Inputs
%    x    1-d signal at fine scale
%    sqmf symmetric filter
%  Outputs
%    y    1-d signal at coarse scale
%
%  See Also
%    DownDyadLo_PBS, UpDyadHi_PBS, UpDyadLo_PBS, FWT_PBS, symm_iconv
%
d = symm_iconv( MirrorSymmFilt(qmf),lshift(x));
n = length(d);
d = d(1:2:(n-1));

%
% Copyright (c) 1995. David L. Donoho
%




%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%

function y = MirrorSymmFilt(x)
% MirrorSymmFilt -- apply (-1)^t modulation to symmetric filter
%  Usage
%    h = MirrorSymmFilt(l)
%  Inputs
%    l   symmetric filter
%  Outputs
%    h   symmetric filter with DC frequency content shifted
%        to Nyquist frequency
%
%  Description
%    h(t) = (-1)^t  * x(t),  -k <= t <= k ; length(x)=2k+1
%
%  See Also
%    DownDyadHi_PBS
%
k = (length(x)-1)/2;
y = ( (-1).^((-k):k) ) .*x;

%
% Copyright (c) 1993. Iain M. Johnstone
%


%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%

function y = symm_iconv(sf,x)
% symm_iconv -- Symmetric Convolution Tool for Two-Scale Transform
%  Usage
%    y = iconv(sf,x)
%  Inputs
%    sf  symmetric filter
%    x   1-d signal
%  Output
%    y    filtered result
%
%  Description
%    Filtering by periodic convolution of x with sf
%
%  See Also
%    symm_aconv, UpDyadHi_PBS, UpDyadLo_PBS, DownDyadHi_PBS, DownDyadLo_PBS
%
n = length(x);
p = length(sf);
if p <= n,
    xpadded = [x((n+1-p):n) x];
else
    z = zeros(1,p);
    for i=1:p,
        imod = 1 + rem(p*n -p + i-1,n);
        z(i) = x(imod);
    end
    xpadded = [z x];
end
ypadded = filter(sf,1,xpadded);
y = ypadded((p+1):(n+p));

shift = (p+1)/2;
shift = 1 + rem(shift-1, n);
y = [y(shift:n) y(1:(shift-1))];



%
% Copyright (c) 1995. Shaobing Chen and David L. Donoho
%



%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%

function x = IWT_SBS(wc,L,qmf,dqmf)
% iwt_po -- Inverse Wavelet Transform (symmetric extension, biorthogonal, symmetric)
%  Usage
%    x = IWT_SBS(wc,L,qmf,dqmf)
%  Inputs
%    wc     1-d wavelet transform: length(wc)= 2^J.
%    L      Coarsest scale (2^(-L) = scale of V_0); L << J;
%    qmf     quadrature mirror filter
%    dqmf    dual quadrature mirror filter (symmetric, dual of qmf)
%  Outputs
%    x      1-d signal reconstructed from wc
%  Description
%    Suppose wc = FWT_SBS(x,L,qmf,dqmf) where qmf and dqmf are orthonormal
%    quad. mirror filters made by MakeBioFilter.  Then x can be reconstructed
%    by
%      x = IWT_SBS(wc,L,qmf,dqmf)
%  See Also:
%    FWT_SBS, MakeBioFilter
%

wcoef = ShapeAsRow(wc);
[n,J] = dyadlength(wcoef);

dp = dyadpartition(n);

x = wcoef(1:dp(L+1));

for j=L:J-1,
    dyadj = (dp(j+1)+1):dp(j+2);
    x = UpDyad_SBS(x, wcoef(dyadj), qmf, dqmf);
end
x = ShapeLike(x,wc);

%
% Copyright (c) 1996. Thomas P.Y. Yu
%

%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%

function x = UpDyad_SBS(beta,alpha,qmf,dqmf)
% UpDyad_SBS -- Symmetric Upsampling operator
%  Usage
%    x = UpDyad_SBS(beta,alpha,qmf,dqmf)
%  Inputs
%    beta  coarse coefficients
%    alpha fine coefficients
%    qmf   quadrature mirror filter
%    dqmf  dual quadrature mirror filter
%  Outputs
%    x     1-d signal at fine scale
%  See Also
%    DownDyad_SBS, IWT_SBS
%

% oddf = (rem(length(qmf),2)==1);
oddf = ~(qmf(1)==0 & qmf(length(qmf))~=0);
oddx = (length(beta) ~= length(alpha));

L = length(beta)+length(alpha);

if oddf,
    if oddx,
        ebeta = extend(beta,1,1);
        ealpha = extend(alpha,2,2);
    else
        ebeta = extend(beta,2,1);
        ealpha = extend(alpha,1,2);
    end
else
    if oddx,
        ebeta = extend(beta,1,2);
        ealpha = [alpha 0 -reverse(alpha)];
    else
        ebeta = extend(beta,2,2);
        ealpha = [alpha -reverse(alpha)];
    end
end

coarse = UpDyadLo_PBS(ebeta,dqmf);
fine = UpDyadHi_PBS(ealpha,qmf);
x = coarse + fine;
x = x(1:L);

%
% Copyright (c) 1996. Thomas P.Y. Yu
%


%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%

function y = UpDyadLo_PBS(x,qmf)
% UpDyadLo_PBS -- Lo-Pass Upsampling operator; periodized
%  Usage
%    u = UpDyadLo_PBS(d,sf)
%  Inputs
%    d    1-d signal at coarser scale
%    sf   symmetric filter
%  Outputs
%    u    1-d signal at finer scale
%
%  See Also
%    DownDyadLo_PBS , DownDyadHi_PBS , UpDyadHi_PBS, IWT_PBS, symm_iconv
%
y =  symm_iconv(qmf, UpSample(x,2) );

%
% Copyright (c) 1995. David L. Donoho
%




%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%

function y = UpDyadHi_PBS(x,qmf)
% UpDyadHi_PBS -- Hi-Pass Upsampling operator; periodized
%  Usage
%    u = UpDyadHi_PBS(d,f)
%  Inputs
%    d    1-d signal at coarser scale
%    sf   symmetric filter
%  Outputs
%    u    1-d signal at finer scale
%
%  See Also
%    DownDyadLo_PBS, DownDyadHi_PBS, UpDyadLo_PBS, IWT_PBS, symm_aconv
%

y = symm_aconv( MirrorSymmFilt(qmf), rshift( UpSample(x,2) ) );

%
% Copyright (c) 1995. David L. Donoho
%





%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%


function wc = FWT2_SBS(x,L,qmf,dqmf)
% FWT2_SBS -- 2-dimensional wavelet transform
%              (symmetric extension, bi-orthogonal)
%  Usage
%    wc = FWT2_SBS(x,L,qmf,dqmf)
%  Inputs
%    x     2-d image (n by n array, n arbitrary)
%    L     coarsest level
%    qmf   low-pass quadrature mirror filter
%    dqmf  high-pass dual quadrature mirror filter
%  Output
%    wc    2-d wavelet transform
%  Description
%    A two-dimensional Wavelet Transform is computed for the
%    matrix x. To reconstruct, use IWT2_SBS.
%  See Also
%    IWT2_SBS
%

[m,J] = dyadlength(x(:,1));
[n,K] = dyadlength(x(1,:));
wc = x;
mc = m;
nc = n;

J = min([J,K]);

for jscal=J-1:-1:L,

    if rem(mc,2)==0,
        top = (mc/2+1):mc;
        bot = 1:(mc/2);
    else
        top = ((mc+1)/2+1):mc;
        bot = 1:((mc+1)/2);
    end
    if rem(nc,2)==0,
        right = (nc/2+1):nc;
        left = 1:(nc/2);
    else
        right = ((nc+1)/2+1):nc;
        left = 1:((nc+1)/2);
    end

    for ix=1:mc,
        row = wc(ix,1:nc);
        [beta,alpha] = DownDyad_SBS(row,qmf,dqmf);
        wc(ix,left) = beta;
        wc(ix,right) = alpha;
    end
    for iy=1:nc,
        column = wc(1:mc,iy)';
        [beta,alpha] = DownDyad_SBS(column,qmf,dqmf);
        wc(bot,iy) = beta';
        wc(top,iy) = alpha';
    end
    mc = bot(length(bot));
    nc = left(length(left));
end







%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%

function x = IWT2_SBS(wc,L,qmf,dqmf)
% IWT2_SBS -- Inverse 2d Wavelet Transform
%            (symmetric extention, bi-orthogonal)
%  Usage
%    x = IWT2_SBS(wc,L,qmf,dqmf)
%  Inputs
%      wc    2-d wavelet transform [n by n array, n arbitrary]
%      L     coarse level
%      qmf   low-pass quadrature mirror filter
%      dqmf  high-pas dual quadrature mirror filter
%  Outputs
%      x     2-d signal reconstructed from wc
%  Description
%      If wc is the result of a forward 2d wavelet transform, with
%           wc = FWT2_SBS(x,L,qmf,dqmf)
%      then x = IWT2_SBS(wc,L,qmf,dqmf) reconstructs x exactly if qmf is a nice
%      quadrature mirror filter, e.g. one made by MakeBioFilter
%  See Also:
%    FWT2_SBS, MakeBioFilter
%

[m,J] = dyadlength(wc(:,1));
[n,K] = dyadlength(wc(1,:));
% assume m==n, J==K

x = wc;

dpm = dyadpartition(m);

for jscal=L:J-1,
    bot = 1:dpm(jscal+1);
    top = (dpm(jscal+1)+1):dpm(jscal+2);
    all = [bot top];

    nc = length(all);

    for iy=1:nc,
        x(all,iy) =  UpDyad_SBS(x(bot,iy)', x(top,iy)', qmf, dqmf)';
    end
    for ix=1:nc,
        x(ix,all) = UpDyad_SBS(x(ix,bot), x(ix,top), qmf, dqmf);
    end
end

%
% Copyright (c) 1996. Thomas P.Y. Yu
%

%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%

function wc = FWT2_PO(x,L,qmf)
% FWT2_PO -- 2-d MRA wavelet transform (periodized, orthogonal)
%  Usage
%    wc = FWT2_PO(x,L,qmf)
%  Inputs
%    x     2-d image (n by n array, n dyadic)
%    L     coarse level
%    qmf   quadrature mirror filter
%  Outputs
%    wc    2-d wavelet transform
%
%  Description
%    A two-dimensional Wavelet Transform is computed for the
%    array x.  To reconstruct, use IWT2_PO.
%
%  See Also
%    IWT2_PO, MakeONFilter
%
[n,J] = quadlength(x);
wc = x;
nc = n;
for jscal=J-1:-1:L,
    top = (nc/2+1):nc; bot = 1:(nc/2);
    for ix=1:nc,
        row = wc(ix,1:nc);
        wc(ix,bot) = DownDyadLo(row,qmf);
        wc(ix,top) = DownDyadHi(row,qmf);
    end
    for iy=1:nc,
        row = wc(1:nc,iy)';
        wc(top,iy) = DownDyadHi(row,qmf)';
        wc(bot,iy) = DownDyadLo(row,qmf)';
    end
    nc = nc/2;
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

function x = IWT2_PO(wc,L,qmf)
% IWT2_PO -- Inverse 2-d MRA wavelet transform (periodized, orthogonal)
%  Usage
%    x = IWT2_PO(wc,L,qmf)
%  Inputs
%    wc    2-d wavelet transform [n by n array, n dyadic]
%    L     coarse level
%    qmf   quadrature mirror filter
%  Outputs
%    x     2-d signal reconstructed from wc
%
%  Description
%    If wc is the result of a forward 2d wavelet transform, with
%    wc = FWT2_PO(x,L,qmf), then x = IWT2_PO(wc,L,qmf) reconstructs x
%    exactly if qmf is a nice qmf, e.g. one made by MakeONFilter.
%
%  See Also
%    FWT2_PO, MakeONFilter
%
[n,J] = quadlength(wc);
x = wc;
nc = 2^(L+1);
for jscal=L:J-1,
    top = (nc/2+1):nc; bot = 1:(nc/2); all = 1:nc;
    for iy=1:nc,
        x(all,iy) =  UpDyadLo(x(bot,iy)',qmf)'  ...
            + UpDyadHi(x(top,iy)',qmf)';
    end
    for ix=1:nc,
        x(ix,all) = UpDyadLo(x(ix,bot),qmf)  ...
            + UpDyadHi(x(ix,top),qmf);
    end
    nc = 2*nc;
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


%
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%



function wp = FWT_TI(x,L,qmf)
% FWT_TI -- translation invariant forward wavelet transform
%  Usage
%    TIWT = FWT_TI(x,L,qmf) 
%  Inputs
%    x        array of dyadic length n=2^J
%    L        degree of coarsest scale
%    qmf      orthonormal quadrature mirror filter 
%  Outputs
%    TIWT     stationary wavelet transform table
%             formally same data structure as packet table
%
%  See Also
%    IWT_TI
%

	[n,J] = dyadlength(x);
	D = J-L;
	wp = zeros(n,D+1);
	x = ShapeAsRow(x);
%
	wp(:,1) = x';
	for d=0:(D-1),
		for b=0:(2^d-1),
		   s = wp(packet(d,b,n),1)';
		   hsr = DownDyadHi(s,qmf);
		   hsl = DownDyadHi(rshift(s),qmf);
		   lsr = DownDyadLo(s,qmf);
		   lsl = DownDyadLo(rshift(s),qmf);
		   wp(packet(d+1,2*b  ,n),d+2) = hsr';
		   wp(packet(d+1,2*b+1,n),d+2) = hsl';
		   wp(packet(d+1,2*b  ,n),1  ) = lsr';
		   wp(packet(d+1,2*b+1,n),1  ) = lsl';		   
		 end
	end

%
% Copyright (c) 1994. David L. Donoho
% 
    
    
%   
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%   
    

function tiwt = FWT2_TI(x,L,qmf)
% FWT_TI -- 2-D translation invariant forward wavelet transform
%  Usage
%    TIWT = FWT2_TI(x,L,qmf) 
%  Inputs
%    x        2-d image (n by n real array, n dyadic)
%    L        degree of coarsest scale
%    qmf      orthonormal quadrature mirror filter 
%  Outputs
%    TIWT     translation-invariant wavelet transform table, (3*(J-L)+1)*n by n
%
%  See Also
%    IWT2_TI, IWT2_TIMedian
%

	[n,J] = quadlength(x);
	D = J-L;

	tiwt = zeros((3*D+1)*n, n);
	lastx = (3*D*n+1):(3*D*n+n); lasty = 1:n;
	tiwt(lastx,lasty) = x;
%
	for d=0:(D-1),
	  l = J-d-1; ns = 2^(J-d);
	  for b1=0:(2^d-1), for b2=0:(2^d-1),
	      s = tiwt(3*D*n+packet(d,b1,n), packet(d,b2,n));

	      wc00 = FWT2_PO(s,l,qmf);
	      wc01 = FWT2_PO(CircularShift(s,0,1),l,qmf);
	      wc10 = FWT2_PO(CircularShift(s,1,0),l,qmf);
	      wc11 = FWT2_PO(CircularShift(s,1,1),l,qmf);

	      index10 = packet(d+1,2*b1,n); index20 = packet(d+1,2*b2,n);
	      index11 = packet(d+1,2*b1+1,n); index21 = packet(d+1,2*b2+1,n);
	      % horizontal stuff
	      tiwt(3*d*n + index10 , index20) = wc00(1:(ns/2),(ns/2+1):ns);
	      tiwt(3*d*n + index11,  index20) = wc01(1:(ns/2),(ns/2+1):ns);
	      tiwt(3*d*n + index10 , index21) = wc10(1:(ns/2),(ns/2+1):ns);
	      tiwt(3*d*n + index11 , index21) = wc11(1:(ns/2),(ns/2+1):ns);
	      % vertical stuff
	      tiwt((3*d+1)*n + index10 , index20) = wc00((ns/2+1):ns,1:(ns/2));
	      tiwt((3*d+1)*n + index11,  index20) = wc01((ns/2+1):ns,1:(ns/2));
	      tiwt((3*d+1)*n + index10 , index21) = wc10((ns/2+1):ns,1:(ns/2));
	      tiwt((3*d+1)*n + index11 , index21) = wc11((ns/2+1):ns,1:(ns/2));
	      % diagonal stuff
	      tiwt((3*d+2)*n + index10 , index20) = wc00((ns/2+1):ns,(ns/2+1):ns);
	      tiwt((3*d+2)*n + index11,  index20) = wc01((ns/2+1):ns,(ns/2+1):ns);
	      tiwt((3*d+2)*n + index10 , index21) = wc10((ns/2+1):ns,(ns/2+1):ns);
	      tiwt((3*d+2)*n + index11 , index21) = wc11((ns/2+1):ns,(ns/2+1):ns);
	      % low freq stuff
	      tiwt(3*D*n + index10 , index20) = wc00(1:(ns/2),1:(ns/2));
	      tiwt(3*D*n + index11,  index20) = wc01(1:(ns/2),1:(ns/2));
	      tiwt(3*D*n + index10 , index21) = wc10(1:(ns/2),1:(ns/2));
	      tiwt(3*D*n + index11 , index21) = wc11(1:(ns/2),1:(ns/2));
	    end, end
	end

% 
% Copyright (c) 1995. David L. Donoho and Thomas P.Y. Yu
% 
    
    
%   
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%   
    
function p=packet(d,b,n)
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
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%   
    

function x = IWT_TI(pkt,qmf)
% IWT_TI -- Invert translation invariant wavelet transform
%  Usage
%    x = IWT_TI(TIWT,qmf)
%  Inputs
%    TIWT     translation-invariant wavelet transform table
%    qmf      quadrature mirror filter
%  Outputs
%    x        1-d signal reconstructed from translation-invariant
%             transform TIWT
%
%  See Also
%    FWT_TI
%
	[n,D1] = size(pkt); 
	D = D1-1;
	J = log2(n);
	L = J-D;
%
	wp = pkt;
%
	sig = wp(:,1)'; 
	for d= D-1:-1:0,  
		for b=0:(2^d-1)
			hsr = wp(packet(d+1,2*b  ,n),d+2)';
		    hsl = wp(packet(d+1,2*b+1,n),d+2)';
		    lsr = sig(packet(d+1,2*b  ,n) );
		    lsl = sig(packet(d+1,2*b+1,n) );		   
			loterm = (UpDyadLo(lsr,qmf) + lshift(UpDyadLo(lsl,qmf)))/2;
			hiterm = (UpDyadHi(hsr,qmf) + lshift(UpDyadHi(hsl,qmf)))/2;
			sig(packet(d,b,n)) = loterm+hiterm;
		end
	end
	x = sig;

%
% Copyright (c) 1994. David L. Donoho
% 
    
    
%   
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%   
    


function x = IWT2_TI(tiwt,L,qmf)
% IWT2_TI -- Invert 2-d translation invariant wavelet transform
%  Usage
%    x = IWT2_TI(TIWT,qmf)
%  Inputs
%    TIWT     translation-invariant wavelet transform table, (3*(J-L)+1)*n by n
%    L        degree of coarsest scale
%    qmf      quadrature mirror filter
%  Outputs
%    x        2-d image reconstructed from translation-invariant transform TIWT
%
%  See Also
%    FWT2_TI, IWT2_TIMedian
%
	[D1,n] = size(tiwt);
	J = log2(n);
	D = J-L;
%
	lastx = (3*D*n+1):(3*D*n+n); lasty = 1:n;
	x = tiwt(lastx,lasty);

	for d=(D-1):-1:0,
	  l = J-d-1; ns = 2^(J-d);
	  for b1=0:(2^d-1), for b2=0:(2^d-1),
	      index10 = packet(d+1,2*b1,n); index20 = packet(d+1,2*b2,n);
	      index11 = packet(d+1,2*b1+1,n); index21 = packet(d+1,2*b2+1,n);

	      wc00 = [x(index10,index20) tiwt(3*d*n+index10,index20) ; ...
		tiwt((3*d+1)*n+index10,index20) tiwt((3*d+2)*n+index10,index20)];
	      wc01 = [x(index11,index20) tiwt(3*d*n+index11,index20) ; ...
		tiwt((3*d+1)*n+index11,index20) tiwt((3*d+2)*n+index11,index20)];
	      wc10 = [x(index10,index21) tiwt(3*d*n+index10,index21) ; ...
		tiwt((3*d+1)*n+index10,index21) tiwt((3*d+2)*n+index10,index21)];
	      wc11 = [x(index11,index21) tiwt(3*d*n+index11,index21) ; ...
		tiwt((3*d+1)*n+index11,index21) tiwt((3*d+2)*n+index11,index21)];

	      x(packet(d,b1,n), packet(d,b2,n)) = ( IWT2_PO(wc00,l,qmf) + ....
		CircularShift(IWT2_PO(wc01,l,qmf),0,-1) + ...
		CircularShift(IWT2_PO(wc10,l,qmf),-1,0) + ...
		CircularShift(IWT2_PO(wc11,l,qmf),-1,-1) ) / 4;
	  end, end
	end

% 
% Copyright (c) 1995. David L. Donoho and Thomas P.Y. Yu
%     
    
%   
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%   
    


function result = CircularShift(matrix, colshift, rowshift)
% CIRCULARSHIFT: Circular shifting of a matrix/image, i.e., pixels that get
% shifted off one side of the image are put back on the other side.
%
% result = circularShift(matrix, colshift, rowshift)
% 
% EPS, DJH '96

lastrow = size(matrix, 1);
lastcol = size(matrix, 2);

result = matrix;

% Shift the cols
if (colshift>0)
  result = [result(:,[lastcol-colshift+1:lastcol]) ...
	    result(:,[1:lastcol-colshift])];
else
  colshift = -colshift;
  result = [result(:,[colshift+1:lastcol]) ...
	    result(:,[1:colshift])];
end

% Shift the rows
if (rowshift>0)
  result = [result([lastrow-rowshift+1:lastrow],:) ; ...
	    result([1:lastrow-rowshift],:)];
else
  rowshift = -rowshift;
  result = [result([rowshift+1:lastrow],:) ; ...
	    result([1:rowshift],:)];
end

function x = reverse(x)

x = x(end:-1:1);


function StatWT = TI2Stat(TIWT)
% TI2Stat -- Convert Translation-Invariant Transform to Stationary Wavelet Transform
%  Usage
%    StatWT = TI2Stat(TIWT)
%  Inputs
%    TIWT     translation invariant table from FWT_TI
%  Outputs
%    StatWT   stationary wavelet transform table table as FWT_Stat
%
%  See Also
%    Stat2TI, FWT_TI, FWT_Stat
%
	StatWT = TIWT; 
	[n,D1] = size(StatWT); 
	D = D1-1;
	J = log2(n);
	L = J-D;
%
	index = 1;
	
	for d=1:D,
		nb = 2^d;
		nk = n/nb;
		
		index = [ (index+nb/2); index];
		index = index(:)';
		
		for b= 0:(nb-1),
			StatWT(d*n + (index(b+1):nb:n)) = TIWT(d*n + packet(d,b,n));
		end
	end
	
	for b= 0:(nb-1),
		StatWT((index(b+1):nb:n)) = TIWT(packet(d,b,n));
	end

%
% Copyright (c) 1994. Shaobing Chen
%

function TIWT = Stat2TI(StatWT)
% Stat2TI -- Convert Stationary Wavelet Transform to Translation-Invariant Transform 
%  Usage
%    TIWT = Stat2TI(StatWT)
%  Inputs
%    StatWT  stationary wavelet transform table as FWT_Stat
%  Outputs
%    TIWT    translation-invariant transform table as FWT_TI
%
%  See Also
%    Stat2TI, FWT_TI, FWT_Stat
%

	TIWT = StatWT; 
	[n,D1] = size(StatWT); 
	D = D1-1;
	index = 1;
	
	for d=1:D,
		nb = 2^d;
		nk = n/nb;
		
		index = [ (index+nb/2); index];
		index = index(:)';
		
		for b= 0:(nb-1),
			TIWT(d*n + packet(d,b,n))  = StatWT(d*n + (index(b+1):nb:n));
		end
	end
	
	for b= 0:(nb-1),
		 TIWT(packet(d,b,n)) = StatWT((index(b+1):nb:n));
	end

%
% Copyright (c) 1994. Shaobing Chen
%

	
	
    
     
 
%
%  Part of Wavelab Version 850
%  Built Tue Jan  3 13:20:40 EST 2006
%  This is Copyrighted Material
%  For Copying permissions see COPYING.m
%  Comments? e-mail wavelab@stat.stanford.edu 
