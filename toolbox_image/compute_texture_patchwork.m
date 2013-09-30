function [M,Id] = compute_texture_patchwork(H,n, options)

% compute_texture_patchwork - mix several textures
%
%   M = compute_texture_patchwork(H,n,options);
%
%   H can be a cell array of textures or a (n,n,3) matrix.
%   Handles up to 5 textures.
%
%   To generate a random partition, set options.patchwork_mode='random'.
%
%   Copyright (c) 2006 Gabriel Peyre

if not(iscell(H))
    H1 = H; H = {};
    for i=1:size(H1,3)
        H{i} = H1(:,:,i);
    end
end

options.null = 0;
mode = getoptions(options, 'patchwork_mode', 'deterministic');

if nargin<2
    n = size(H{1},1);
end



p = length(H); % number of textures

if strcmp(mode, 'random')
        
    sigma = getoptions(options, 'patchwork_sigma', 60);
    J = floor(log2(p));
    A = ones(n);
    for j=0:J-1
        B = A;
        for k=1:2^j
            I = find(B==k);
            U = perform_blurring(randn(n),sigma,options);
            s = median(U(I));
            I1 = find( (B==k) & (U>s) );
            I2 = find( (B==k) & (U<=s) );
            A(I1) = 2*k-1;
            A(I2) = 2*k;
        end
    end
    
    M = zeros(n);
    Id = zeros(n);
    for i=1:p
        B = H{i};
        if size(B,1)<n
            B = perform_image_extension(B,n);
        end
        B = B(1:n,1:n);
        M(A==i) = B(A==i);
        Id(A==i) = i;
    end
    return;

end

m = min(size(H{1},1),n);

M = zeros(n);
M(1:m,1:m) = H{1}(1:m,1:m);
Id = ones(n);

if p==2
    M(end/2+1:end,:) = H{2}(end-n/2+1:end,:);
    Id(end/2+1:end,:) = 2;
    return;
end

if p>1
    M(1:end/2,end/2+1:end) = H{2}(1:1:n/2,1:n/2);
    Id(1:end/2,end/2+1:end) = 2;
end
if p>2
    M(end/2+1:end,1:end/2) = H{3}(1:1:n/2,1:n/2);
    Id(end/2+1:end,1:end/2) = 3;
end
if p>3
    M(end/2+1:end,end/2+1:end) = H{4}(1:1:n/2,1:n/2);
    Id(end/2+1:end,end/2+1:end) = 4;
end
if p>4
    r = 1;
    A = H{5}(1:1:n/2,1:n/2);
    x = linspace(-1,1,n/2)';
    [Y,X] = meshgrid(x,x);
    J = find( X.^2 + Y.^2 <= r^2 );
    x  = [ones(n/4,1)*Inf; x; ones(n/4,1)*Inf];
    [Y,X] = meshgrid(x,x);
    I = find( X.^2 + Y.^2 <= r^2 );
    M(I) = A(J);
    Id(I) = 5;
end