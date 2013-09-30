function M = perform_varying_blurring(M, Sigma, options)

% perform_varying_blurring - perform a spacially varying gaussian blurring
%
%   M = perform_varying_blurring(M, Sigma, options);
%
%   Sigma(i,j) is the width of blurring (in pixels) around pixel (i,j).
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
if size(M,3)>1
    for i=1:size(M,3)
        M(:,:,i) = perform_varying_blurring(M(:,:,i), Sigma, options);
    end
    return;
end

n = size(M,1);
if isfield(options, 'nb_filters')
    m = nb_filters;
else
    m = 15;
end

if std(Sigma(:))<1e-8
    % special case: in fact constant blurring
    m = 1;
end

sigma = linspace( min(Sigma(:)), max(Sigma(:)), m );
Mh = zeros(n,n,m);
for i=1:m
    Mh(:,:,i) = perform_blurring(M, sigma(i), options);
end

% combine together the filters
[tmp,I] = min( abs( repmat(Sigma,[1 1 m]) - repmat( reshape(sigma, [1 1 m]), [n n 1]) ), [], 3 );
J = (1:n^2) + (I(:)'-1)*n^2; J = reshape(J, n,n);
M = Mh(J);