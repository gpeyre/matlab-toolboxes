function M1 = perform_conditional_histogram_matching(M_src,M_src_cond, M_tgt,M_tgt_cond, nb_bins, nb_bins_cond)

% compute_conditional_histogram - compute conditional histograms
%
% M1 = perform_conditional_histogram_matching(M_src,M_src_cond, M_tgt,M_tgt_cond, , nb_bins, nb_bins_cond)
%
%   Match the histograms of M_src (conditionned to M_src_cond) to the histograms 
%   of M_tgt conditionned to M_tgt_cond.
%
%   Copyright (c) 2005 Gabriel Peyré


if nargin<5
    nb_bins = 21;
end
if nargin<6
    nb_bins_cond = 21;
end

M1 = M_src;

% compute targte histograms
[H,HX,HX_cond] = compute_conditional_histogram(M_tgt,M_tgt_cond, nb_bins, nb_bins_cond);

bins_size_cond = HX_cond(2)-HX_cond(1);

for i=1:length(HX_cond)
    v = HX_cond(i);
    % locate the coefficient to transform
    I = find( abs(M_src_cond-v)<bins_size_cond/2 );
    x = M_src(I);
    % retrieve the histograms and bins centers
    X = HX{i};
    N = H(i,:);
    if ~isempty(X) & ~isempty(x)
        % perform the matching
        warning off;
        x = histoMatch(x, N, X);
        warning on;
        % assign the values
        M1(I) = x;
    end
end