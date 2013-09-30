function M_src = perform_spaciowav_matching(M_src, M_tgt, nb_bins, options)

% perform_spaciowav_matching - match both the spacial and wavelet histogram
%
% M_src = perform_spaciowav_matching(M_src, M_tgt, nb_bins, options);
%
%   Copyright (c) 2005 Gabriel Peyré

options.null = 0;
if nargin<3
    nb_bins = 100;
end

%% for cell array
if iscell(M_src)
    for i=1:length(M_src)
        M_src{i} = perform_spaciowav_matching(M_src{i}, M_tgt{i}, nb_bins, options);
    end
    return;
end

if ~isfield(options, 'decomp_type')
    options.decomp_type = 'tri';
end
if ~isfield(options, 'wavelet_vm')
    options.wavelet_vm = 3;
end
if ~isfield(options, 'wavelet_type')
    options.wavelet_type = 'biorthogonal';
end
if isfield(options, 'callback_transform')
    callback_transform = options.callback_transform;
else
    callback_transform = @perform_atrou_transform;
end

n = size(M_src,1);
Jmax = log2(n)-1;
Jmin = max(Jmax-5, 2);

% spacial matching
M_src = perform_histogram_matching(M_src, M_tgt, nb_bins);

% perform forward transform
M_src_wav = feval(callback_transform,M_src,Jmin,options);
M_tgt_wav = feval(callback_transform,M_tgt,Jmin,options);
p = length(M_src_wav);

% wavelet domain matching
for i=1:p
    M_src_wav{i} = perform_histogram_matching(M_src_wav{i}, M_tgt_wav{i}, nb_bins);
end

% un-transform
M_src = feval(callback_transform,M_src_wav,Jmin,options);

% spacial matching
M_src = perform_histogram_matching(M_src, M_tgt, nb_bins);