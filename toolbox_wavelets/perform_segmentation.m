function [B,err] = perform_segmentation(E,options)

% perform_segmentation - perform image segmentation
%
%   B = perform_segmentation(E,options);
%
%   E is an (n,n,k) set of k dimensional features vectors (one per pixel in
%   the image).
%
%   options.segmentation_method can be 
%       'simple': E(:,:,k) should be minimum in area where texture i is
%           present. The algorithm uses a simple filtering.
%       'chan-vese': E(:,:,k) should be minimum in area where texture i is
%          present. The algorithm uses a Chan-Vese active contour. 
%       'kmeans': E(i,j,:) is a feature vector around pixel (i,j).
%           the algorithm compute feature center for each texture using
%           k-means and then perform nearest-neighbors classification.
%
%   Copyright (c) 2007 Gabriel Peyre

n = size(E,1);

ntextures = size(E,3);

segmentation_method = getoptions(options, 'segmentation_method', 'simple');
oracle = getoptions(options, 'oracle', []);
lambda = getoptions(options, 'lambda', 1.2);

% for testing various filtering sizes
if isempty(oracle)
    nmu = 1;
    mu = getoptions(options, 'mu', 4);
    mu_list = mu;
else
    nmu = getoptions(options, 'nmu', 18);
    mumax = getoptions(options, 'mumax', 10);
    mu_list = linspace( 0.1,mumax,nmu );
end

% best dictionary
switch segmentation_method
    case 'simple'
        
        q = size(E,3);
        q1 = min(20,q); % PCA dimension
        err_list = []; B1 = {};
        for i=1:nmu
            progressbar(i,nmu);
            EE = perform_blurring(E,mu_list(i));
            [tmp,B1{i}] = min(EE,[],3);
            if not(isempty(oracle))
                err_list(i) = 1 - sum( B1{i}(:)~=oracle(:) )/n^2;
            else
                err_list(i) = 0;
            end
        end
        [err,i] = max(err_list);
        B = B1{i};
        
    case 'kmeans'
        
        %% perform clustering via Kmeans
        ntextures = getoptions(options, 'ntextures',[],1);
        q = size(E,3);
        q1 = min(20,q); % PCA dimension
        err_list = []; B1 = {};
        for i=1:nmu
            progressbar(i,nmu);
            % apply non-linearity + smoothing
            E1 = perform_blurring(E,mu_list(i));
            E1 = reshape(E1,[n^2 q])';
            % add spacial relationship
            % rescale the features
            E1 = E1 - repmat(min(E1,[],2),[1 n^2]);
            E1 = E1 ./ repmat( median(E1,2),[1 n^2]);
            % dimension reduction
            [Y,E1] = pca(E1,q1);
            % do k-means on a subset
            options.nb_iter = 50;
            options.kmeans_code = 3;
            [B1{i},seeds,rms] = perform_kmeans(E1,ntextures,options);
            B1{i} = reshape(B1{i},n,n);
            [B1{i},err_list(end+1)] = fit_segmentation(B1{i},oracle);            
        end
        [err,i] = max(err_list);
        B = B1{i};
            
    case 'chan-vese'

        % active contours
        En = zeros(n,n,ntextures);
        if ntextures==2
            En = E(:,:,1)-E(:,:,2);
        elseif ntextures==4
            v = {[4 2] [3 1] [4 3] [2 1] };
            for k=1:4
                En(:,:,k) = E(:,:,v{k}(1))-E(:,:,v{k}(2)); a =  En(:,:,k);
                % En(:,:,k) = En(:,:,k)-median(a(:)); a =  En(:,:,k);
                % En(:,:,k) = lambda * En(:,:,k)/std( a(:) );
            end
        else
            error('Works only for 2/4 textures');
        end
        % normalize
        options.E = lambda * En/sqrt( mean(En(:).^2) );
        if not(isfield(options, 'Tmax'))
            options.Tmax = 180;
        end
        options.redistance_freq = 10;
        options.dt = 0.5;
        options.display_freq = 10;
        if not(isfield(options, 'nb_svg'))
            options.nb_svg = 0;
        end
        options.solver = 'grad';
        namec = 'small-disks';  % initial shape
        Dist0 = compute_levelset_shape(namec, n);
        if ntextures==4 % shift this function for the other phase
            sel = mod( (1:n) + 8,n )+1;
            Dist0 = cat(3,Dist0,Dist0(sel,sel));
        end
        fprintf('Active contours:');
        Dist = perform_active_contour(Dist0, 'chan-vese-user', options);
        if ntextures==2
            B = (Dist>0)+1;
        else
            B = zeros(n); D1 = Dist(:,:,1); D2 = Dist(:,:,2);
            B(D1<0 & D2<0) = 1;
            B(D1<0 & D2>=0) = 2;
            B(D1>=0 & D2<0) = 3;
            B(D1>=0 & D2>=0) = 4;
        end
end

err = 0;
if not(isempty(oracle))
    [B,err] = fit_segmentation(B,oracle);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [B,err] = fit_segmentation(B,B0)

if isempty(B0)
    return;
end

n = size(B,1);
% check for optimal permutation to fit B0
ntextures = max(B(:));
P = perms(1:ntextures);
err1 = [];
for k=1:size(P,1)
    pp = P(k,:)';
    err1(k) = sum( P(k,B(:))'==B0(:) )/n^2;
end
[err,k] = max(err1);
B = reshape( P(k,B), size(B) );