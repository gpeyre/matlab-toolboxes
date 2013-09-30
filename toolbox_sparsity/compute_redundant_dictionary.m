function [D,info] = compute_redundant_dictionary(name,n,options)

% compute_redundant_dictionary - compute several redundant dictionaries (matrices)
%
%   [D,info] = compute_redundant_dictionary(name,n,options);
%
%   n is the dimension
%   options.q controls the redundancy
%   info is a struct continaining information about the dictionary.
%
%   name can be 'wavelets' (with redundant translation)
%   'cosine' (with redundant frequencies)
%   'gabor'
%
%   for Gabor, you can select options.gaborloc which controls the ratio of
%   frequency/time localization. gaborloc==2 means that the gabor atoms are
%   2x larger localization in frequencies.
%
%   Copyright (c) 2008 Gabriel Peyre

q = getoptions(options, 'q', 4);

options.null = 0;
switch name
    case 'wavelets'
        [D,info.sel] = compute_wav_dictionary(n,q,options);
    case 'cosine'
        D = compute_cosine_dictionary(n,q,options);
        info.sel = 1:size(D,2);
    case 'gabor'
        [D,info] = compute_gabor_dictionary(n,q,options);
        info.sel = 1:size(D,2);
    otherwise
        error('Unknown dictionary.');
end


%%%
function [D,info] = compute_gabor_dictionary(n,q,options)


%% compute the freq/time localization factor
gaborloc = getoptions(options, 'gaborloc', 1);
gaborcpx = getoptions(options, 'gaborcpx', 0);
[sigma,sigmaf] = compute_gabor_sigma(n, gaborloc);
info.sigma = sigma;
info.sigmaf = sigmaf;

%% build the windowing functions
t = [0:n/2 -n/2+1:-1];
g = exp(-t.^2 / (2*sigma^2) );

sf = sqrt( sigmaf/sigma *n/q ); % spacing in time
sx = sqrt( sigma /sigmaf*n/q ); % spacing in frequency

x = ( 0:sx:n-1 );
if gaborcpx
    f = ( 0:sf:n-1 );
else
    f = ( 0:sf:n/2-1 );
end
t = 0:n-1;

nx = length(x);
nf = length(f);
p = nx*nf; % number of atoms

%% compute gabor atoms
[T,F,X] = ndgrid( t, f, x );
TX = T-X; % shifted
if gaborcpx
    D = exp( -TX.^2 / (2*sigma^2) ) .* exp( 2i*pi/n * TX .* F );
else
    D = exp( -TX.^2 / (2*sigma^2) ) .* cos( 2*pi/n * TX .* F );
end
D = reshape(D, [n p]);
D = D ./ repmat( sqrt(sum(abs(D).^2)), [n 1] );

% frequency index
F = F(1,:,:); info.F = F(:);
% spacial index
X = X(1,:,:); info.X = X(:);
info.nx = nx; info.nf = nf;



%% select the sigma so that the localization freq/time
function [sigma,sigmaf] = compute_gabor_sigma(n, gaborloc)

% is gaborloc
ns = 400;
slist = linspace(1,n/8,ns);
t = [0:n/2 -n/2+1:-1];
[S,T] = meshgrid( slist,t );
% generate gaussians
G = exp(-T.^2 ./ (2*S.^2) );
G = G ./ repmat( sqrt(sum(G.^2)), [n 1] );
% compute FFT
GF = real(fft(G));
GF = GF ./ repmat( sqrt(sum(GF.^2)), [n 1] );
% correlation
C = G'*GF;
% C(i,j) : correlation space_i / freq_j
[tmp,I] = max(C, [], 2);
slistf = slist(I);
r = slistf./slist;
[tmp,I] = min(abs(r-gaborloc)); I = I(1);
if I==1 || I==ns
    warning('Problem in detecting Gabor scale');
end
sigma = slist(I);
sigmaf = slistf(I);

%%%
function D = compute_cosine_dictionary(n,q,options)

s = 0:1/q:n/2; s(end) = []; 
[F,X] = meshgrid(s,0:n-1);
D = cos( 2*pi/n * F.*X);  

%%%
function [D,selj] = compute_wav_dictionary(n,q,options)

Jmax = log2(n);
Jmin = 2;
options.wavelet_type = 'biorthogonal_swapped';
options.wavelet_vm = 4;

D = [];
sel = [n/2+1:n];
selj = {};
for j=Jmax-2:-1:Jmin
    y = zeros(n,1);
    y(sel(end/2)) = 1;
    psi = perform_wavelet_transform(y, Jmin-1, -1, options);
    sel = sel(1:end/2); sel = sel-length(sel);
    % translate
    qj = max(1/q*n/2^j,1);
    
    [Y,X] = meshgrid(1:qj:n,1:n);
    D = [D psi(mod(X-Y,n)+1)];
    selj{j} = size(D,2)-n/qj+1:size(D,2);
end

