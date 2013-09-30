function T = compute_tensor_field_random(n,options)

% compute_tensor_field_random - create 2D TF
%
%   T = compute_tensor_field_random(n,options);
%
%   Copyright (c) 2008 Gabriel Peyre

options.null = 0;
sigma_flow = getoptions(options, 'sigma_tensor', 50*n/256);
niter_tensor = getoptions(options, 'niter_tensor', 8);
verb = getoptions(options, 'verb', 1);

T = randn(n,n,3);
en = linspace(.5,1,n^2);
an = linspace(.1,1,n^2);
or = linspace(-pi/2,pi/2,n^2);
for i=1:niter_tensor
    if verb
        progressbar(i,niter_tensor);
    end
    T = perform_blurring(T,sigma_flow);
    U = perform_tensor_mapping(T,+1);
%    clf;
%   for k=1:3
%        subplot(1,3,k); a = U(:,:,k);
%        hist(a(:),100); axis tight;
%    end
%    drawnow;
    U(:,:,1) = perform_histogram_equalization(U(:,:,1), en);
    U(:,:,2) = perform_histogram_equalization(U(:,:,2), an);
    U(:,:,3) = perform_histogram_equalization(U(:,:,3), or);
    T = perform_tensor_mapping(U,-1);
end