function [y,quincunx_filters] = perform_quincunx_wavelet_transform(x,Jmin,dir,options)

% perform_quincunx_wavelet_transform - compute quincunx transform
%
% Forward transform
%   [MW,options.quincunx_filters] = perform_quincunx_wavelet_transform(M,Jmin,+1,options);
% Backward transform
%   M = perform_quincunx_wavelet_transform(MW,Jmin,-1,options);
%
% This is a simple wrapper to the code of
%   Dimitri Van De Ville
%   Biomedical Imaging Group, BIO-E/STI
%   Swiss Federal Institute of Technology Lausanne
%   CH-1015 Lausanne EPFL, Switzerland
%
% For more information, see
%   Isotropic Polyharmonic B-Splines: Scaling Functions and Wavelets
%   D. Van De Ville, T. Blu, M. Unser
%   IEEE Transactions on Image Processing, vol. 14, no. 11, pp. 1798-1813, November 2005.
%
% options.type is the type of transform:
%    - McClellan: O/B/D,
%    - Polyharmonic Rabut: Po/Pb/Pd
%    - Polyharmonic isotropic: PO/PB/PD
% options.gamma is the degree of the wavelet transform

options.null = 0;


% Setup parameters of wavelet transform
% -------------------------------------
% 1. Type of wavelet transform
%    - McClellan: O/B/D,
%    - Polyharmonic Rabut: Po/Pb/Pd
%    - Polyharmonic isotropic: PO/PB/PD
if isfield(options, 'type')
    type = options.type;
else
    type='PO';
end

% 2. Degree of the wavelet transform
%    (order gamma polyharmonic=degree alpha+2)
if isfield(options, 'gamma')
    gamma = options.gamma;
else
    gamma = 4;
end
alpha=gamma-2;

% 3. Number of iterations
n = size(x,1);
Jmax = log2(n)-1;
J = 2*(Jmax-Jmin+1);


% Precalculate filters
% --------------------
if isfield(options, 'quincunx_filters')
    quincunx_filters = options.quincunx_filters;
    FA = quincunx_filters{1};
    FS = quincunx_filters{2};
else
    [FA,FS]=FFTquincunxfilter2D([size(x)],alpha,type);
    quincunx_filters{1} = FA;
    quincunx_filters{2} = FS;
end

% Perform wavelet analysis
% ------------------------
if dir==1
    y = FFTquin2D_analysis(x,J,FA); y=real(y);
else
    y = FFTquin2D_synthesis(x,J,FS);
end


function A0=autocorr2d(H)

% AUTOCORR2D
%       Iterative frequency domain calculation of the 2D autocorrelation
%       function. A=autocorr2d(H) computes the frequency response of
%       the autocorrelation filter A(exp(2*i*Pi*nu)) corresponding to the
%       scaling function with refinement filter H.
%       Please note that the 2D grid of the refinement filter (nu),
%       corresponds to a uniformly sampled grid, twice as coarse
%       (in each dimension) as the given filter H.
%
% Please treat this source code as confidential!
% Its content is part of work in preparation to be submitted:
%
% T. Blu, D. Van De Ville, M. Unser, ''Numerical methods for the computation
% of wavelet correlation sequences,'' SIAM Numerical Analysis.
%
% Biomedical Imaging Group, EPFL, Lausanne, Switzerland.


len1=size(H,1)/2;
len2=size(H,2)/2;

% stop criterion
crit=1e-4;
improvement=1.0;

% initial "guess"
A0=ones(len1,len2);

% calculation loop

while improvement>crit,

    % sinc interpolation
    if 1,
        Af=fftshift(ifftn(A0));
        Af(1,:)=Af(1,:)/2;
        Af(:,1)=Af(:,1)/2;
        Af(len1+1,1:len2)=Af(1,len2:-1:1);
        Af(1:len1,len2+1)=Af(len1:-1:1,1);
        Ai=zeros(2*len1,2*len2);
        Ai(len1-len1/2+1:len1+len1/2+1,len2-len2/2+1:len2+len2/2+1)=Af;
        Ai=fftn(ifftshift(Ai));
    else
        Ai=interp2(A0,1,'nearest'); Ai(:,end+1)=Ai(:,end); Ai(end+1,:)=Ai(end,:);
    end;

    % recursion
    A1=Ai(1:len1,1:len2).*H(1:len1,1:len2)+Ai(len1+1:2*len1,1:len2).*H(len1+1:2*len1,1:len2)+Ai(1:len1,len2+1:2*len2).*H(1:len1,len2+1:2*len2)+Ai(len1+1:2*len1,len2+1:2*len2).*H(len1+1:2*len1,len2+1:2*len2);

    improvement=mean(abs(A1(:)-A0(:)));

    A0=A1;

end;


function y2=FFTquin2D_analysis(x,J,FA);

% FFTQUIN2D_ANALYSIS Perform 2D quincunx wavelet analysis
% 	y2 = FFTquin2D_analysis(x,J,FA)
%
% 	Input:
%       x  = input data
%       J  = number of iterations
%       FA = analysis filters
%
%       Output:
%       y2 = result (folded)
%
% Dimitri Van De Ville
%      Biomedical Imaging Group, BIO-E/STI
%      Swiss Federal Institute of Technology Lausanne
%      CH-1015 Lausanne EPFL, Switzerland

% Dimensions of data
[m,n]=size(x);

% Dimensions of remaining folded lowpass subband
cm=m; cn=n;

% Check dimensions of input data
for iter=1:2,
    tmp=size(x,iter)/2^floor((J+2-iter)/2);
    if round(tmp)~=tmp,
        disp(' ')
        disp('The size of the input signal must be a sufficient power of two!')
        disp(' ')
        y2=[];
        return;
    end;
end;

% Analysis filters:
%  H1,   G1   : filters for iterations with mod(j,2)=1
%  H1D,  G1D  : filters for iterations with mod(j,2)=0
% (premultiply analysis filters by the downsampling factor; i.e., det(D)=2)
H1=FA(:,:,1)/2;
G1=FA(:,:,2)/2;
H1D=FA(:,:,3)/2;
G1D=FA(:,:,4)/2;

% Compute the FFT of the input data
X=fftn(x); clear x;

% Initialize cube for folded wavelet coefficients
y2=zeros(cm,cn);

for j=1:J,

    mj=mod(j,2);

    % Select filter
    switch mj,
        case 1,
            H=H1;
            G=G1;
        case 0,
            H=H1D;
            G=G1D;
    end;

    % Filtering
    Y1=H(1:size(X,1),1:size(X,2)).*X;
    Y2=G(1:size(X,1),1:size(X,2)).*X;

    % Downsampling & upsampling
    switch mj,
        case 1,
            if j~=J,
                % This operates like: Y1=Y1+fftshift(Y1); Y1=Y1(1:cm,1:cn/2);
                Y1(1:cm/2,1:cn/2)=Y1(1:cm/2,1:cn/2)+Y1(cm/2+1:cm,cn/2+1:cn);
                Y1(cm/2+1:cm,1:cn/2)=Y1(cm/2+1:cm,1:cn/2)+Y1(1:cm/2,cn/2+1:cn);
                Y1=Y1(1:cm,1:cn/2);
            else
                Y1=Y1+fftshift(Y1);
            end;
            Y2=Y2 + fftshift(Y2);
            %! y2(1:cm,cn/2+1:cn)=fold2D(real(ifftn(Y2)),mj,m,n);
            y2(1:cm,cn/2+1:cn)=fold2D((ifftn(Y2)),mj,m,n);
            cn=cn/2;
        case 0,
            m2=m/2;n2=n/2;
            Y1=Y1(1:m2,1:n2)+Y1(m2+1:2*m2,1:n2);
            Y2=Y2(1:m2,1:n2)+Y2(m2+1:2*m2,1:n2);
            y2(cm/2+1:cm,1:cn)=real(ifftn(Y2));
            cm=cm/2;

            % Preparing filters for next three iterations
            lm=1:2:m; ln=1:2:n;
            H1=H1(lm,ln);
            G1=G1(lm,ln);
            H1D=H1D(lm,ln);
            G1D=G1D(lm,ln);
            m=m2; n=n2;
    end;

    X=Y1;
end;

% Insert folded lowpass subband
y2(1:cm,1:cn)=fold2D(real(ifftn(Y1)),mj,m,n);
return;


% ----------------------------------------------------------------------
% Fold 2D subband
% ----------------------------------------------------------------------
function f=fold2D(u,j,m,n)
switch j,
    case 0,
        f=u;
    case 1,
        f=u(:,1:2:n)+u(:,2:2:n);
end;


function x=FFTquin2D_synthesis(y2,J,FS);

% FFTQUIN2D_SYNTHESIS Perform 2D quincunx wavelet synthesis
% 	x = FFTquin2D_synthesis(y2,J,FS)
%
% 	Input:
%       y2 = data (folded)
%       J = number of interations
%       FS = synthesis filters
%
%       Output:
%       x = result (unfolded)
%
% Dimitri Van De Ville
%      Biomedical Imaging Group, BIO-E/STI
%      Swiss Federal Institute of Technology Lausanne
%      CH-1015 Lausanne EPFL, Switzerland

% Dimensions of input data
[m,n]=size(y2);

% Synthesis filters:
%  H1,   G1   : filters for iterations with mod(j,2)=1
%  H1D,  G1D  : filters for iterations with mod(j,2)=0
H1=FS(:,:,1);
G1=FS(:,:,2);
H1D=FS(:,:,3);
G1D=FS(:,:,4);

% How many "double-iterations"?
l=2^((J-mod(J,2))/2);
M=m/l; N=n/l;

% Dimensions of full data
[M0,N0,P0]=size(H1);

% Dimensions of folded lowpass subband
cm=M0; cn=N0;
cn=cn/2^floor((J+1)/2);
cm=cm/2^floor((J)/2);

% Extract folded lowpass subband and unfold
y1=y2(1:cm,1:cn);
y1=unfold2D(y1,mod(J,2),M,N);
Y1=fftn(y1); clear y1;

j=J;
while j>0;
    mj=mod(j,2);

    % unfold highpass subband and select filters
    switch mj,
        case 1,
            f=unfold2D(y2(1:M,N/2+1:N),mj,M,N);
            H=H1(1:l:M0,1:l:N0);
            G=G1(1:l:M0,1:l:N0);
        case 0,
            f=unfold2D(y2(M+1:2*M,1:N),mj,M,N);
            l=l/2;
            Y1=[Y1 Y1;Y1 Y1];
            H=H1D(1:l:M0,1:l:N0);
            G=G1D(1:l:M0,1:l:N0);
    end;

    Y2=fftn(f); clear f;
    if mj==0,
        Y2=[Y2 Y2;Y2 Y2];
        M=2*M; N=2*N;
    end;

    % filter data
    Y1=Y1.*H + Y2.*G;

    j=j-1;
end;
x=real(ifftn(Y1));
return;

% ----------------------------------------------------------------------
% Unfold 2D subband
% ----------------------------------------------------------------------
function u=unfold2D(y,j,M,N)
switch j,
    case 0,
        u=y;
    case 1,
        u=zeros(M,N);
        u(1:2:M,1:2:N)=y(1:2:M,:);
        u(2:2:M,2:2:N)=y(2:2:M,:);
end;
return;


function [FA,FS]=FFTquincunxfilter2D(dim,alpha,type);

% FFTQUINCUNXFILTER2D Supplies filters for the 2D quincunx transform.
% 	[FFTanalysisfilters,FFTsynthesisfilters]=FFTquincunxfilter2D(dim,alpha,type)
% 	computes the frequency response of low- and high-pass filters.
%
% 	Input:
% 	dim = dimensions of the input signal
% 	alpha = degree of the filters, any real number >=0
%   type =
%     class I: McClellan-based filters
%       o = orthogonal; b = discrete linear B-spline (synthesis); d = dual
%     class II: Polyharmonic B-spline wavelets (Rabut-style localization)
%       Po = orthogonal; Pb = B-spline; Pd = dual
%     class III: Polyharmonic B-spline wavelets (isotropic-style
%                localization)
%       PO = orthogonal {=Po}; PB = B-spline; PD = dual
%
% Dimitri Van De Ville
%      Biomedical Imaging Group, BIO-E/STI
%      Swiss Federal Institute of Technology Lausanne
%      CH-1015 Lausanne EPFL, Switzerland

m=dim(1);
n=dim(2);

gamma=alpha+2;

% Setup downsampling matrix
D=[1 1; 1 -1]';

% Coordinates in the Fourier domain
[xo,yo]=ndgrid(2*pi*([1:m]-1)/m,2*pi*([1:n]-1)/n);

for iter=0:1,

    if iter==0, % First filter: no downsampling
        x=xo;
        y=yo;
    end;
    if iter>0,  % Coordinates after downsampling
        x=D(1,1)*xo+D(1,2)*yo;
        y=D(2,1)*xo+D(2,2)*yo;
        % D=D*D;    % prepare filter for next iteration
    end;

    if lower(type(1)) == 'p', % Isotropic polyharmonic splines

        if iter==0,
            [ac,acD,loc,B] = FFTquincunxpolyfilter(x,y,gamma,type(2));
        else
            %ac = fftshift(interp2(xo,yo,ac,mod(x,2*pi),mod(y,2*pi),'*nearest'),1);
            %acD= fftshift(interp2(xo,yo,acD,mod(x,2*pi),mod(y,2*pi),'*nearest'));
            %loc= fftshift(interp2(xo,yo,loc,mod(x,2*pi),mod(y,2*pi),'*nearest'),1);
            %B  = fftshift(interp2(xo,yo,B,mod(x,2*pi),mod(y,2*pi),'*nearest'),1);
            ac = interp2(xo,yo,ac0,mod(x,2*pi),mod(y,2*pi),'*nearest');
            acD= interp2(xo,yo,acD0,mod(x,2*pi),mod(y,2*pi),'*nearest');
            loc= interp2(xo,yo,loc0,mod(x,2*pi),mod(y,2*pi),'*nearest');
            B  = interp2(xo,yo,B0,mod(x,2*pi),mod(y,2*pi),'*nearest');
            %[ac,acD,loc,B] = FFTquincunxpolyfilter(x,y,gamma,type(2));
        end;
        ortho = acD./ac;

    else

        [ortho,B] = FFTquincunxmcclellan(x,y,alpha);

    end;

    % Lowpass filter
    [H,H1] = FFTquincunxfilter2D_lowpass(x,y,type,B,ortho);

    % Reversed filter
    B0 = B;
    if iter==0,
        B=fftshift(B); ortho=fftshift(ortho);
    else
        B=fftshift(B,1); ortho=fftshift(ortho,1);
    end;

    if lower(type(1)) == 'p',
        ac0=ac; acD0=acD; loc0=loc;
        if iter==0,
            ac=fftshift(ac);
            acD=fftshift(acD);
            loc=fftshift(loc);
        else
            ac=fftshift(ac,1);
            acD=fftshift(acD);
            loc=fftshift(loc,1);
        end;
    else
        ac=0; ac0=0; acD=0; acD0=0; loc=0; loc0=0;
    end;

    [Hd,H1d] = FFTquincunxfilter2D_highpass(x-pi,y-pi,type,B,ortho,ac,ac0,acD,acD0,loc,loc0);

    % Highpass filters: modulation
    G  = exp(i*x).*H1d;
    G1 = exp(-i*x).*Hd;

    % Fill up array with analysis and synthesis filters
    FA(:,:,iter*2+1) = H1;
    FA(:,:,iter*2+2) = G1;
    FS(:,:,iter*2+1) = H;
    FS(:,:,iter*2+2) = G;
end;



function [H,H1]=FFTquincunxfilter2D_highpass(x,y,type,B,ortho,ac,ac0,acD,acD0,loc,loc0);

% FFTQUINCUNXFILTER2D_HIGHPASS Supplies highpass filters for 2D quincunx transform.
% 	[Hsynthesis,Hanalysis]=FFTquincunxfilter2D_highpass(x,y,alpha,type)
% 	computes the frequency response of highpass filters.
%
% 	Input:
% 	x,y = coordinates in the frequency domain
% 	alpha = degree of the filters, any real number >=0
%       type = type of filter (ortho, dual, bspline)
%
% Dimitri Van De Ville
%      Biomedical Imaging Group, BIO-E/STI
%      Swiss Federal Institute of Technology Lausanne
%      CH-1015 Lausanne EPFL, Switzerland


if lower(type(1))=='p', % Polyharmonic splines

    type = [ type(2:length(type)) ];

    switch lower(type(1)),

        % Orthogonal filters
        case 'o',
            H = B*sqrt(2) ./ sqrt(ortho);
            H1  = conj(H);

            % Dual filters (B-spline at the analysis side)
        case 'd',
            H = sqrt(2)*conj(B).*ac;
            H1  = sqrt(2)*B./acD;

            % B-spline filters (B-spline at the synthesis side)
        case 'b',
            H1 = sqrt(2)*conj(B).*ac;
            H  = sqrt(2)*B./acD;

        case 'i',
            H = sqrt(2)*loc0./ac0;
            H1 = sqrt(2)*(B.^2).*ac./acD.*ac0./loc0;
            H1(find(isinf(H1)))=0;

    end;

else  % Simple fractional iterate with McClellan transform
    [H,H1]=FFTquincunxfilter2D_lowpass(x,y,type,B,ortho);
end;



function [H,H1]=FFTquincunxfilter2D_lowpass(x,y,type,B,ortho);

% FFTQUINCUNXFILTER2D_LOWPASS Supplies lowpass filters for 2D quincunx transform.
% 	[Hsynthesis,Hanalysis]=FFTquincunxfilter2D_lowpass(x,y,alpha,type)
% 	computes the frequency response of lowpass filters.
%
% 	Input:
% 	x,y = coordinates in the frequency domain
% 	alpha = degree of the filters, any real number >=0
%       type = type of filter (ortho, dual, bspline)
%
% Dimitri Van De Ville
%      Biomedical Imaging Group, BIO-E/STI
%      Swiss Federal Institute of Technology Lausanne
%      CH-1015 Lausanne EPFL, Switzerland

H=B;

if lower(type(1)) == 'p',
    type = [ type(2:length(type)) ];
end;

switch lower(type(1)),

    % Orthogonal filters
    case 'o',
        H = H*sqrt(2) ./ sqrt(ortho);
        H1  = conj(H);

        % Dual filters (B-spline at the analysis side)
    case 'd',
        H1 = sqrt(2)*conj(H);
        H = conj(H1./(ortho));

        % B-spline filters (B-spline at the synthesis side)
    case 'b',
        H = sqrt(2)*H;
        H1 = conj(H./(ortho));

        % Isotropic wavelet filter
    case 'i',
        H1 = sqrt(2)*conj(H);
        H = conj(H1./(ortho));

end;

function [ortho,H] = FFTquincunxmcclellan(x,y,alpha)

% FFTQUINCUNXMCCLELLAN Supplies orthonormalizing part of the 2D McClellan
%       filters for 2D quincunx transform.
%
% 	v=FFTquincunmcclellan(x,y,alpha)
%
% 	Input:
% 	x,y = coordinates in the frequency domain
% 	alpha = degree of the filters, any real number >=1
%
% Dimitri Van De Ville
%      Biomedical Imaging Group, BIO-E/STI
%      Swiss Federal Institute of Technology Lausanne
%      CH-1015 Lausanne EPFL, Switzerland


% McClellan transform filters
H  = (2*ones(size(x)) + (cos(x) + cos(y))).^((alpha+1)/2);
H  = H/4^((alpha+1)/2);
Hc = (2*ones(size(x)) - (cos(x) + cos(y))).^((alpha+1)/2);
Hc = Hc/4^((alpha+1)/2);

% Orthonormalizing denominator (l_2 dual)
ortho = H.*conj(H) + Hc.*conj(Hc);

function [ac,acD,loc,B] = FFTquincunxpolyfilter(x,y,gamma,type)

% FFTQUINCUNXPOLYFILTER Supplies orthonormalizing part of the 2D polyharmonic
%       spline filters for 2D quincunx transform.
%
% 	v=FFTquincunxpolyfilter(x,y,gamma,type)
%
% 	Input:
% 	x,y = coordinates in the frequency domain
% 	gamma = order of the filters, any real number >=0
%
% Dimitri Van De Ville
%      Biomedical Imaging Group, BIO-E/STI
%      Swiss Federal Institute of Technology Lausanne
%      CH-1015 Lausanne EPFL, Switzerland


% Get dimensions
[m n]=size(x);
m=2*m; n=2*n;

% Construct fine grid [0,2pi[ x [0,2pi[
[xo,yo]=ndgrid(2*pi*([1:m]-1)/m,2*pi*([1:n]-1)/n);

warning off MATLAB:divideByZero;

% Refinement filter for [2 0; 0 2]
if upper(type(1)) == type(1),
    loc = 8/3*(sin(xo/2).^2+sin(yo/2).^2)+2/3*(sin((xo+yo)/2).^2+sin((xo-yo)/2).^2);
    H  = 2^(-gamma)*((8/3*(sin(xo).^2+sin(yo).^2)+2/3*(sin(xo+yo).^2+sin(xo-yo).^2))./loc).^(gamma/2);
else
    loc = sin(xo/2).^2+sin(yo/2).^2;
    H  = 2^(-gamma)*((sin(xo).^2+sin(yo).^2)./loc).^(gamma/2);
end;
H(find(isnan(H))) = 1;

% Autocorrelation function
ac = autocorr2d(H.^2);

% Construct fine grid [0,2pi] x [0,2pi]
[xo2,yo2]=ndgrid(2*pi*([1:m/2+1]-1)/(m/2),2*pi*([1:n/2+1]-1)/(n/2));

% Extend autocorrelation function
ac(m/2+1,:)=ac(1,:);
ac(:,n/2+1)=ac(:,1);

% Compute autocorrelation function on subsampled grid; D=[1 1; 1 -1]
x2=mod(xo2+yo2,2*pi);
y2=mod(xo2-yo2,2*pi);

% Compute values on grid (x,y)
%------------------------------
% Autocorrelation filter
acD= interp2(xo2,yo2,ac,mod(x+y,2*pi),mod(x-y,2*pi),'*nearest');
ac = interp2(xo2,yo2,ac,mod(x,2*pi),mod(y,2*pi),'*nearest');

% Localization operator
loc = interp2(xo,yo,loc,mod(x,2*pi),mod(y,2*pi),'*nearest');
loc = loc.^(gamma/2);

% Scaling filter for quincunx lattice
if upper(type(1)) == type(1),
    B = 2^(-gamma/2)*((8/3*(sin((x+y)/2).^2+sin((x-y)/2).^2)+2/3*(sin(x).^2+sin(y).^2))./(8/3*(sin(x/2).^2+sin(y/2).^2)+2/3*(sin((x+y)/2).^2+sin((x-y)/2).^2))).^(gamma/2);
else
    B = 2^(-gamma/2)*((sin((x+y)/2).^2+sin((x-y)/2).^2)./(sin(x/2).^2+sin(y/2).^2)).^(gamma/2);
end;
B(find(isnan(B))) = 1;

function [v,fv] = poly(alpha,type)

% POLY - compute the polyharmonic B-spline
%
% 	Input:
% 	alpha = degree of the polyharmonic spline
%       type = type of B-spline (b,d,o,B,D,O)
%
% Dimitri Van De Ville
%      Biomedical Imaging Group, BIO-E/STI
%      Swiss Federal Institute of Technology Lausanne
%      CH-1015 Lausanne EPFL, Switzerland


warning off MATLAB:divideByZero;

N=16; % spatial support: [-N/4:N/4[
Z=16; % spatial zoom factor
%N=64;
%Z=2;

step=2*pi/N;
[x,y]=meshgrid(-Z*pi:step:Z*pi-step);
w1=x; w2=y;

% compute frequency response
gamma=alpha+2;

x = x/2; y = y/2;
loc = sin(x).^2+sin(y).^2;

% alternative localization
if upper(type(1)) == type(1),
    loc = ( 8/3*loc + 2/3*(sin(x+y).^2+sin(x-y).^2) ) / 4;
end;

pow=(x.^2+y.^2);

loc=loc.^(gamma/2);
pow=pow.^(gamma/2);

fv  = ( loc ./ pow );
fv(find(isnan(fv))) = 1;

% autocorrelation function
[xo,yo]=meshgrid(0:step/2:2*pi-step/2);

% refinement filter for [2 0; 0 2]
if upper(type(1)) == type(1),
    H  = 2^(-gamma)*((8/3*(sin(xo).^2+sin(yo).^2)+2/3*(sin(xo+yo).^2+sin(xo-yo).^2))./(8/3*(sin(xo/2).^2+sin(yo/2).^2)+2/3*(sin((xo+yo)/2).^2+sin((xo-yo)/2).^2))).^(gamma/2);
else
    H  = 2^(-gamma)*((sin(xo).^2+sin(yo).^2)./(sin(xo/2).^2+sin(yo/2).^2)).^(gamma/2);
end;
H(find(isnan(H))) = 1;

% autocorrelation function
ac0 = autocorr2d(H.^2);

sx = size(ac0,1); sy = size(ac0,2);
ac = zeros(size(fv));
for iterx=1:Z,
    for itery=1:Z,
        bx=(iterx-1)*sx+1;
        by=(itery-1)*sy+1;
        ac(bx:bx+sx-1,by:by+sy-1)=ac0;
    end;
end;
ac = fftshift(ac);

if lower(type(1)) == 'o',
    fv = fv ./ sqrt(ac);
end;

if lower(type(1)) == 'd',
    fv = fv ./ ac;
end;

% compute spatial version
fv = ifftshift(fv);              % shift to [0:2pi[
v = real((ifft2(fv)))*Z^2;       % compensate for zooming

% only keep the half [-N/4:N/4[ instead of [-N/2:N/2[ (to avoid aliasing)
[x,y]=meshgrid(-N/4:1/Z:N/4-1/Z);
v=fftshift(v);
v=v((N*Z)/4:(N*Z)*3/4-1,(N*Z)/4:(N*Z)*3/4-1);
v=ifftshift(v);

% optional: normalize energy (for movies)
if length(type)>1 & type(2) == 'n',
    s=sum(v(:).^2)/(Z*Z);
    v=v*sqrt(1/s);
end;

% surface plot with lighting
surfl(x,y,fftshift(v)); shading flat; view(-26,44);
xlabel('x_1');
ylabel('x_2');

function [v,fv] = polyw(alpha,type)

% POLYW - polyharmonic B-spline wavelet
%
% 	Input:
%
%       type
%         lowercase - Rabut style
%         uppercase - isotropic brand
%
% 	alpha = degree of the polyharmonic spline
%
% Dimitri Van De Ville
%      Biomedical Imaging Group, BIO-E/STI
%      Swiss Federal Institute of Technology Lausanne
%      CH-1015 Lausanne EPFL, Switzerland


gamma=alpha+2;

warning off MATLAB:divideByZero;

N=16;
Z=16;
step=2*pi/N;
[w1,w2]=meshgrid(-Z*pi:step:Z*pi-step);

% downsampled grid
w1D=(w1+w2)/2;
w2D=(w1-w2)/2;

% scaling function
w1D  = w1D/2; w2D = w2D/2;

loc  = 4*sin(w1D).^2+4*sin(w2D).^2;

if upper(type(1)) == type(1),
    loc = loc - 8/3*sin(w1D).^2.*sin(w2D).^2;
end;

pow = (2*w1D).^2+(2*w2D).^2;

loc=loc.^(gamma/2);
pow=pow.^(gamma/2);

fv1  = loc ./  pow;
fv1(find(isnan(fv1))) = 1;
w1D  = 2*w1D; w2D = 2*w2D;

% autocorrelation function
[w1o,w2o]=meshgrid(0:step/4:2*pi-step/4);

% refinement filter for [2 0; 0 2]
if upper(type(1)) == type(1),
    H  = 2^(-gamma)*((8/3*(sin(w1o).^2+sin(w2o).^2)+2/3*(sin(w1o+w2o).^2+sin(w1o-w2o).^2))./(8/3*(sin(w1o/2).^2+sin(w2o/2).^2)+2/3*(sin((w1o+w2o)/2).^2+sin((w1o-w2o)/2).^2))).^(gamma/2);
else
    H  = 2^(-gamma)*((sin(w1o).^2+sin(w2o).^2)./(sin(w1o/2).^2+sin(w2o/2).^2)).^(gamma/2);
end;
H(find(isnan(H))) = 1;

% autocorrelation function
ac0 = autocorr2d(H.^2);
ac0(size(ac0,1)+1,:)=ac0(1,:);
ac0(:,size(ac0,2)+1)=ac0(:,1);

[w1o,w2o]=meshgrid(0:step/2:2*pi);
ac = interp2(w1o,w2o,ac0,mod(w1D,2*pi),mod(w2D,2*pi),'*linear');

if lower(type(1)) == 'o' || lower(type(1)) == 'j',
    fv1 = fv1 ./ sqrt(ac);
end;
if lower(type(1)) == 'd',
    fv1 = fv1 ./ ac;
end;

% wavelet filter
w1D = -w1D-pi; w2D = -w2D-pi;

acD = interp2(w1o,w2o,ac0,mod(w1D,2*pi),mod(w2D,2*pi),'*linear');

w1D = w1D/2;  w2D = w2D/2;

% Quincunx scaling filter
if upper(type(1)) == type(1),
    t1 = 8/3*(sin(w1D+w2D).^2 + sin(w1D-w2D).^2) + 2/3*(sin(2*w1D).^2 + sin(2*w2D).^2);
    t2 = 8/3*(sin(w1D).^2 + sin(w2D).^2) + 2/3*(sin(w1D+w2D).^2 + sin(w1D-w2D).^2);
else
    t1 = sin(w1D+w2D).^2 + sin(w1D-w2D).^2;
    t2 = sin(w1D).^2 + sin(w2D).^2;
end;

% Dyadic scaling filter
if length(type)>1 & type(2)=='D',
    t1 = sin(2*w1D).^2 + sin(2*w2D).^2;
    t2 = sin(w1D).^2 + sin(w2D).^2;
end;

t=( 0.5 * t1./t2 ).^(gamma/2);
fv2 = conj(t);
fv2(find(isnan(fv2))) = 1;
w1D = 2*w1D; w2D = 2*w2D;
w1D = -w1D-pi; w2D = -w2D-pi;

if lower(type(1)) == 'i',
    fv2 = loc;
end;
if lower(type(1)) == 'j',
    fv2 = loc;
    %  fv2(find(isinf(abs(fv2)))) = 0;
end;

% shift
fv2 = fv2 .* exp(-i*w1D);

acD2 = interp2(w1o,w2o,ac0,mod(w1,2*pi),mod(w2,2*pi),'*linear');

if lower(type(1)) == 'o',
    fv3 = sqrt(acD./acD2);
end;
if lower(type(1)) == 'b',
    fv3 = acD;
end;
if lower(type(1)) == 'd',
    fv3 = 1./acD2;
end;
if lower(type(1)) == 'i',
    fv3 = 1./ac;
end;

% scale relation
fv = 2^(gamma/2-1)*fv1.*fv2.*fv3;
%fv=fv1;
fv = ifftshift(fv);

% compute spatial version
v = ifft2(fv)*Z^2;

[x1,x2]=meshgrid(-N/4:1/Z:N/4-1/Z);
v=fftshift(v);
v=v((N*Z)/4:(N*Z)*3/4-1,(N*Z)/4:(N*Z)*3/4-1);
v=ifftshift(v);

surfl(x1,x2,fftshift(real(v))); shading flat; view(-26,44);
xlabel('x_1');
ylabel('x_2');