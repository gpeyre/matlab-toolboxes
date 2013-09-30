function y = perform_atrou_transform(x,Jmin,options)

% perform_atrou_transform - compute the "a trou" wavelet transform, 
%   i.e. without subsampling.
%
%   w_list = perform_atrou_transform(M,Jmin,options);
%
%   'w_list' is a cell array, w_list{ 3*(j-Jmin)+q }
%   is an imagette of same size as M containing the full transform
%   at scale j and quadrant q.
%   The lowest resolution image is in w_list{3*(Jmax-Jmin)+4} =
%   w_list{end}.
%
%   The ordering is :
%       { H_j,V_j,D_j, H_{j-1},V_{j-1},D_{j-1}, ..., H_1,V_1,D_1, C }
%
%   'options' is an (optional) structure that may contains:
%       options.wavelet_type: kind of wavelet used (see perform_wavelet_transform)
%       options.wavelet_vm: the number of vanishing moments of the wavelet transform.
%
%	When possible, this code tries to use the following fast mex implementation:
%		* Rice Wavelet Toolbox : for Daubechies filters, ie when options.wavelet_type='daubechies'.
%		Only the Windows binaries are included, please refer to 
%			http://www-dsp.rice.edu/software/rwt.shtml
%		to install on other OS.
%
%   You can set decomp_type = 'quad' or  'tri' for the 7-9 transform.
%
%   Copyright (c) 2006 Gabriel Peyr?


% add let it wave a trou transform to the path
% path('./cwpt2/',path);
% add Rice Wavelet Toolbox transform to the path
% path('./rwt/',path);

if nargin<3
    options.null = 0;
end

if isfield(options, 'wavelet_type') 
    wavelet_type = options.wavelet_type;
    wavelet_vm = 4;
else
    wavelet_type = 'biorthogonal';
    wavelet_vm = 3;
end

if isfield(options, 'wavelet_vm')   
    wavelet_vm = options.wavelet_vm;
end

if isfield(options, 'bound')
    bound = options.bound;
else
    bound = 'sym';
end

% first check for LIW interface, since it's the best code
if exist('cwpt2_interface') && ~strcmp(wavelet_type, 'daubechies') &&  ~strcmp(wavelet_type, 'symmlet')
    
    if  strcmp(wavelet_type, 'biorthogonal') || strcmp(wavelet_type, 'biorthogonal_swapped')
        if wavelet_vm~=3
            warning('Works only for 7-9 wavelets.');
        end
        wavelet_types = '7-9';
    elseif strcmp(wavelet_type, 'battle')
        if wavelet_vm~=3
            warning('Works only for cubic spline wavelets.');
        end
        wavelet_types = 'spline';
    else
        error('Unknown transform.');
    end
    
    if isfield(options, 'decomp_type')
        decomp_type = options.decomp_type;
    else
        decomp_type = 'quad';   % either 'quad' or 'tri'
    end
    
    if ~iscell(x)
        dirs = 'forward';
        J = log2(size(x,1)) - Jmin;
        M = x;
    else
        dirs = 'inverse';
        if strcmp(decomp_type, 'tri')
            x = { x{1:end-3}, x{end}, x{end-2:end-1} };
        else
            x = { x{1:end-4}, x{end}, x{end-3:end-1} };
        end
        J = 0;
        M = zeros(size(x{1},1), size(x{1},2), J);
        for i=1:length(x)
            M(:,:,i) = x{i};
        end
    end
    
    y = cwpt2_interface(M, dirs, wavelet_types, decomp_type, J);
    
    if ~iscell(x)
       M = y; y = {};
       for i=1:size(M,3)
           y{i} = M(:,:,i);
       end
       % put low frequency at the end
       if strcmp(decomp_type, 'tri')
           y = { y{1:end-3}, y{end-1:end}, y{end-2} };
       else
           y = { y{1:end-4}, y{end-2:end}, y{end-3} };
       end
    end
    return;
end

% precompute filters
if strcmp(wavelet_type, 'daubechies')
    fname = 'Daubechies';
    if wavelet_vm==1
        fname = 'Haar';
    end 
    qmf = MakeONFilter(fname,wavelet_vm*2);  % in Wavelab, 2nd argument is VM*2 for Daubechies... no comment ...
elseif strcmp(wavelet_type, 'symmlet')
    qmf = MakeONFilter('Symmlet',wavelet_vm);
elseif strcmp(wavelet_type, 'battle')
    qmf = MakeONFilter('Battle',wavelet_vm-1);
elseif strcmp(wavelet_type, 'biorthogonal')
    [qmf,dqmf] = MakeBSFilter( 'CDF', [wavelet_vm,wavelet_vm] );
elseif strcmp(wavelet_type, 'biorthogonal_swapped')
    [dqmf,qmf] = MakeBSFilter( 'CDF', [wavelet_vm,wavelet_vm] );
else
    error('Unknown transform.');
end

g = qmf;                % for phi
gg = g;


if exist('mirdwt') & exist('mrdwt') & strcmp(wavelet_type, 'daubechies')
    %%% USING RWT %%%
    if ~iscell(x)
        n = length(x);
        Jmax = log2(n)-1;
        %%% FORWARD TRANSFORM %%%
        L = Jmax-Jmin+1;
        [yl,yh,L] = mrdwt(x,g,L);
        for j=Jmax:-1:Jmin
            for q=1:3
                s = 3*(Jmax-j)+q-1;
                M = yh(:,s*n+1:(s+1)*n); 
                y{ 3*(j-Jmin)+q } = M;
            end
        end
        y{ 3*(Jmax-Jmin)+4 } = yl;         
    else
        n = length(x{1});
        Jmax = log2(n)-1;
        %%% BACKWARD TRANSFORM %%%
        L = Jmax-Jmin+1;
        if L ~= (length(x)-1)/3
            warning('Jmin is not correct.');
            L = (length(x)-1)/3;
        end
        yl = x{ 3*(Jmax-Jmin)+4 }; 
        yh = zeros( n,3*L*n );
        for j=Jmax:-1:Jmin
            for q=1:3
                s = 3*(Jmax-j)+q-1;
                yh(:,s*n+1:(s+1)*n) = x{ 3*(j-Jmin)+q };
            end
        end
        [y,L] = mirdwt(yl,yh,gg,L);
    end
    return;
end

warning('This code is depreciated.');


% for reconstruction
h = mirror_filter(g);    % for psi
hh = g;
if exist('dqmf')
    gg = dqmf;
    hh = mirror_filter(gg);
end


n = length(x);
Jmax = log2(n)-1;
if iscell(x)
    error('reverse transform is not yet implemented.');    
end
M = x;


Mj = M;     % image at current scale (low pass filtered)

nh = length(h);
ng = length(g);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the transform 
wb = waitbar(0,'Computing a trou transform.');
for j=Jmax:-1:Jmin
    waitbar((Jmax-j)/(Jmax-Jmin),wb);
    
    % 1st put some zeros in between g and h
    nz = 2^(Jmax-j);    % space between coefs
    hj = zeros(nz*(nh-1)+1, 1);
    gj = zeros(nz*(ng-1)+1, 1);
    hj( 1 + (0:nh-1)*nz ) = h;
    gj( 1 + (0:ng-1)*nz ) = g;
    
    %%%% filter on X %%%%
    Mjh = zeros(n,n);
    Mjg = zeros(n,n);
    for i=1:n
        Mjh(:,i) = perform_convolution( Mj(:,i), hj, bound );
        Mjg(:,i) = perform_convolution( Mj(:,i), gj, bound );
    end
    
    %%%% filter on Y %%%%
    Mjhh = zeros(n,n);
    Mjhg = zeros(n,n);
    Mjgh = zeros(n,n);
    Mjgg = zeros(n,n);
    for i=1:n
        Mjhh(i,:) = perform_convolution( Mjh(i,:)', hj, bound )';
        Mjhg(i,:) = perform_convolution( Mjh(i,:)', gj, bound )';
        
        Mjgh(i,:) = perform_convolution( Mjg(i,:)', hj, bound )';
        Mjgg(i,:) = perform_convolution( Mjg(i,:)', gj, bound )';
    end
    
    Mj = Mjgg;
    w_list{ 3*(j-Jmin)+1 } = Mjgh;
    w_list{ 3*(j-Jmin)+2 } = Mjhg;
    w_list{ 3*(j-Jmin)+3 } = Mjhh;
    
end
close(wb);

w_list{ 3*(Jmax-Jmin)+4 } = Mj;

y = w_list;









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





function y = mirror_filter(x)
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