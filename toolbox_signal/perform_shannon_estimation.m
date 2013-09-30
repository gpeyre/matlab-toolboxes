function m = perform_shannon_estimation(M,T)

% evaluate_nbr_bits - evaluate the number of bits to code a signal
%
%   m = perform_shannon_estimation(M);
% OR if you need quantization before
%   m = perform_shannon_estimation(M,T);
% OR if you want to estimate multiple images
%   m = perform_shannon_estimation({M1,M2},T);
%
%   nbr_bit is a lower bound (Shanon entropy bound) on the number 
%       of bits needed to code the wavelet image.
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<2
    T = -1;
end

if iscell(M)
    m = 0;
    for i=1:length(M)
        m = m + perform_shannon_estimation(M{i},T);
    end
    return;
end

m = length(M(:)) * compute_entropy(M,T);