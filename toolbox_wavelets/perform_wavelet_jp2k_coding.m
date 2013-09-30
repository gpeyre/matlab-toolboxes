function [y,nbr_bits] = perform_wavelet_jp2k_coding(x, options)

% perform_wavelet_jp2k_coding - code a wavelet transformed image using JPEG2000.
%
% Coding:
%   [stream,nbr_bits] = perform_wavelet_jp2k_coding(MW, nbr_bits, options);
% decoding:
%   MW = perform_wavelet_jp2k_coding(stream, nbr_bits, options);
%
%   You must provide the original image in options.M_original for coding.
%   You should provide options.Jmin.
%   You can provide options.options.nbr_bits to enforce maximum bit rate.
%
%   Copyright (c) 2005 Gabriel PeyrŽ

% add Jpeg200 toolbox to the path
% path('./jp2k/',path);

options.null = 0;
if isfield(options, 'Jmin')
    Jmin = options.Jmin;
else
    warning('You should provide options.Jmin.')
    Jmin = 2;
end
if size(x,1)>1 && size(x,2)>1
    dir=+1;
else
    dir=-1;
end
if isfield(options, 'bit_depth')
    bit_depth = options.bit_depth;
else
    % assumed coding in 0...255 of the original image
    bit_depth = 8;
end
if isfield(options, 'nbr_bits')
    nbr_bits = options.nbr_bits;
else
    % use 2 bbp as default coding rate
    nbr_bits = size(x,1)*2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dir==1
    if isfield(options, 'M_original')
        M_original = options.M_original;
    else
        M_original = x*0;
    end
    MW1 = perform_jp2_rescaling(x,Jmin,+1);
    warning off;
    M_cropped = int16(M_original);
    warning on;
    % target size in byte
    nbr_bits_header = 0;                                % for large images, it is not important
    target_size = ( nbr_bits_header + nbr_bits ) / 8;   % 8 for bytes instead of bits
    % perform coding
    y = perform_jp2k_encoding(M_cropped, target_size, bit_depth+1, MW1);
else
    [decoded,MW1] = perform_jp2k_encoding(x, 0, 0, 0);
    y = perform_jp2_rescaling(MW1,Jmin,-1);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MW = perform_jp2_rescaling(MW,Jmin,dir)

n = size(MW,1);
Jmax = log2(n)-1;

for j=Jmax:-1:Jmin
    q_min = 1;
    if j==Jmin
        q_min = 0;
    end
    for q=q_min:3
        [selx,sely] = compute_quadrant_selection(j,q);
        MW(selx,sely) = 2^(-(Jmax-j+1)*dir) * MW(selx,sely);
    end
end