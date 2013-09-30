function nbr_bits = perform_shannon_estimation_wavelets(MW, Jmin, T)

% evaluate_nbr_bits_wavelets - evaluate the number of bits to code a wavelet transform
%
%   nbr_bits = perform_shannon_estimation_wavelets(MW, Jmin, T);
%
%   MW is the wavelet transform,
%   Jmin is the minimum scale of the transform,
%   T is the quantification step,
%   nbr_bit is a lower bound (Shanon entropy bound) on the number 
%       of bits needed to code the wavelet image.
%
%   Copyright (c) 2004 Gabriel Peyré

if ndims(MW)==3
    nbr_bits = 0;
    for i=1:size(MW,3)
        nbr_bits = nbr_bits + perform_shannon_estimation_wavelets(MW(:,:,i), Jmin, T);
    end
    return;
end

% create the list of quadrant
MW_list = convert_wavelets2list(MW, Jmin);

nbr_bits = perform_shannon_estimation(MW_list, T);