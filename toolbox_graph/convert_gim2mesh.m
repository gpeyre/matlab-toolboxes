function [vertex,face] = convert_gim2mesh(gim, sub_sample, recenter)

% convert_gim2mesh - convert a geometry image into a mesh 
%
%   [vertex,face] = convert_gim2mesh(gim, sub_sample);
%
%   Copyright (c) 2005 Gabriel Peyré

if nargin<2
    sub_sample=1;
end
if nargin<3
    recenter = 0;
end

if recenter
    for i=1:3
        gim(:,:,i) = gim(:,:,i) - mean(mean(gim(:,:,i)));
    end
end

n = size(gim,1);

if sub_sample>1    
    % first, smooth
    s = 0.01*sub_sample/n;
    k = 5*sub_sample;
    m = ceil(n/100)*2+1;
    h = compute_gaussian_filter([m m],s,[n n]);
    for i=1:3
        % zero padd
        M = gim_symmetric_extension(gim(:,:,i),k);
        M = conv2(M, h, 'same');
        gim(:,:,i) = M(k+1:end-k, k+1:end-k);
    end
    % then sub-sample
    gim2 = zeros(n/sub_sample, n/sub_sample);
    [Y,X] = meshgrid(0:1/(n-1):1, 0:1/(n-1):1);
    [YI,XI] = meshgrid(0:1/(n/sub_sample-1):1, 0:1/(n/sub_sample-1):1);
    for i=1:3
        gim2(:,:,i) = interp2(X',Y',gim(:,:,i),XI',YI');
    end
    gim = gim2;
    
    % zippering to ensure perfect match of boundaries
    gim(1, end:-1:end/2+1,:) = gim(1, 1:end/2,:);
    gim(end, end:-1:end/2+1,:) = gim(end, 1:end/2,:);
    gim(end:-1:end/2+1,1,:) = gim(1:end/2,1,:);
    gim(end:-1:end/2+1,end,:) = gim(1:end/2,end,:);
end

n = size(gim,1);

[vertex,face] = compute_base_mesh('square', log2(n));
for i=1:3
    x = gim(:,:,i);
    vertex(i,:) = x(:);
end


function M1 = gim_symmetric_extension(M, k)

n = size(M,1);

M1 = zeros(n+2*k);
M1(k+1:end-k,k+1:end-k) = M;
% extension
M1(1:k,:) = M1(2*k:-1:k+1,:);
M1(end-k+1:end,:) = M1(end-k:-1:end-2*k+1,:);
M1(:,1:k) = M1(:,2*k:-1:k+1);
M1(:,end-k+1:end) = M1(:,end-k:-1:end-2*k+1);