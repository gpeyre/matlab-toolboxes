function y = perform_tv_correction(x,T)

% perform_tv_correction - perform correction of the image to that it minimizes the TV norm.
%
%   y = perform_tv_correction(x,T);
%
%   Perform correction using thresholding of haar wavelet coefficients on 1
%   scale.
%
%   Copyright (c) 2006 Gabriel Peyre

n = size(x,1);
Jmin = log2(n)-1;

% do haar transform
options.wavelet_type = 'daubechies';
options.wavelet_vm = 1;
y = perform_atrou_transform(x,Jmin,options);

% perform soft thesholding
y = perform_soft_tresholding(y,T);

% do reconstruction
y = perform_atrou_transform(y,Jmin,options);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = perform_soft_tresholding(x,t)

if iscell(x)
    for i=1:length(x)
        y{i} = perform_soft_tresholding(x{i},t);
    end
    return;
end

s = abs(x) - t;
s = (s + abs(s))/2;
y = sign(x).*s;