function [MW1,nbr_bits] = perform_jp2k_degradation(MW,Jmin,nbr_bits,M, bit_depth)

% perform_jp2k_degradation - perform coding/decoding of wavelet coefficients.
%
% MW1 = perform_jp2k_degradation(MW,Jmin,nbr_bits,M, bit_depth);
%
%   MW is the wavelet transform.
%   Jmin is the minimum scale.
%   nbr_bits is the target number of bits of the coded stream.
%   bit_depth is 8 for 0...255 valued images.
%
%   Copyright (c) 2006 Gabriel Peyré

% path('./jp2k/',path);

if nargin<5
    % assumed coding in 0...255 of the original image
    bit_depth = 8;
end

MW1 = perform_jp2_rescaling(MW,Jmin,+1);

warning off;
M_cropped = int16(M);
warning on;

nbr_bits_header = 0;                                % for large images, it is not important
target_size = ( nbr_bits_header + nbr_bits ) / 8;   % 8 for bytes instead of bits

[comp_data] = perform_jp2k_encoding(M_cropped, target_size, bit_depth+1, MW1);
[decoded,MW1] = perform_jp2k_encoding(comp_data, 0, 0, 0);

nbr_bits = length(comp_data)*8;

MW1 = perform_jp2_rescaling(MW1,Jmin,-1);

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