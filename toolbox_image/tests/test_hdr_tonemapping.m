%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test for HDR display
% 
% This is based on the simple global compression of
% luminance proposed in 
%
%   Photographic Tone Reproduction for Digital Images
%   Erik Reinhard, Michael Stark, Peter Shirley, James Ferwerda
%   SIGGRAPH 2002
%
%   Copyright (c) 2005 Gabriel Peyre

rep = '../toolbox_signal_data/hdr/';

name = 'smallOffice.hdr';
name = 'cornellbox.hdr';
name = 'memorial.hdr';


% load image
M = load_hdr([rep name]);
n = size(M,1);

M = crop(M, min(size(M,1),size(M,2)) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transform the image into luminance/chrominance space
H = rgb2ycbcr(M);

% take the luminance (first principal component)
L = H(:,:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP1: set medium luminance

% compute mean value
epsi = 0.01;
m = exp( 1/prod(size(L)) * sum( log(L(:)+epsi)) );
% target medium luminance
a = 0.05;
% L = L*a/m;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP2: compress luminance

% perform global gamma correction
% L1 = perform_gamma_correction(L,'exp'); 
use_gamma = 1;
if use_gamma
    gamma = 0.8;
    L1 = L.^gamma;
else
    % using Reihnard correction curve
    Lwhite = 1.1;   % set Lwhite=1 or Lwhite=max(L(:)) for low dinamyc range
    L1 = L .* (1 + L/Lwhite^2)./(1+L);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% recompose
H1 = H;
H1(:,:,1) = L1;
M1 = ycbcr2rgb(H1);

% clip last percentages
% M1 = histoClip(M1,0.01,0.99);

subplot(1,2,1);
imagesc(clamp(M));
axis image; axis off;
subplot(1,2,2);
imagesc(clamp(M1));
axis image; axis off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gradient domain HDR compression
options.bound = 'sym';
La = log(L);
sigma = 3;
h = compute_gaussian_filter([51 51],sigma/(2*n),[n n]);
La = perform_convolution(La,h,options);
g = grad(La, options);
d = sqrt( sum(g.^2,3) ); d(d<eps)=1;
g = g./repmat(d, [1 1 2]);
% compress gradient
alpha = 0.1*mean(d(:));
beta = 0.85;
% reduction factor
r = (d/alpha).^(beta-1);
d1 = d .* r;
g1 = g.*repmat(d1, [1 1 2]);
% reconstruct using Poisson solver
G = compute_periodic_poisson(div(g1,options),1);
L1 = exp(G);

