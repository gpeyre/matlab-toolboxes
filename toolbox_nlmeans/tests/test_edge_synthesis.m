% test for image synthesis using the manifold model

clear all;
k = 3; % hald window width
n = 100;  % size of synthesized image


%% compute the set of patches
w = 2*k+1;
options.n_delta = 17;
options.n_theta = 16;
options.sigma = 0.05;
options.rescale = 0;
options.manifold_type = 'lines';
options.manifold_type = 'edges';
[Q,delta_list,theta_list] = compute_edge_patches(w, options);
Q = (Q+1)/2; % patches should have values in [0,1]

ma = size(Q,3);
na = 1;
s = 1;  % number of colors

%% PCA projection
options.ndims = min(40,w^2); % number of dimension used for distance computation (PCA dim.reduc.)
H = reshape(Q, [s*w^2 na*ma]);
[options.P,X1,v,options.Psi] = pca(H,options.ndims);
% perform PCA projection
H = H - repmat( options.Psi, [1 na*ma] );
H = options.P'*H;
% reshape matrix
H = reshape(H, [options.ndims ma na]);
H = reshape( shiftdim(H,1), [ma na options.ndims]);
Ma = Q((end+1)/2,(end+1)/2,:); Ma = Ma(:);

options.Ha = H;
options.Ma = Ma;
options.k = k;
    

%% options of NL means
options.T = 0.01;       % width of the gaussian, relative to max(M(:))  (=1 here)
options.max_dist = 100000;  % search width (search all the patches)


% lambda=1 : total synthesis
lambda = 1;


repimg = ['results/' options.manifold_type '-synth/w' num2str(w) '-t' num2str(options.T) '-la' num2str(lambda) '/'];
if not(exist(repimg))
    mkdir(repimg);
end
reptheta = [repimg 'theta-delta/'];
if not(exist(reptheta))
    mkdir(reptheta);
end
name = ['edge'];

% initialization
M = double( rand(n)>0.5 );

use_thresh_decay = 0;

niter = 50;
%% iterations
for i=1:niter
    progressbar(i,niter);

    % set width of the windows
    if use_thresh_decay
        options.T = Tmax - (i-1)/(niter-1)*(Tmax-Tmin);
    end
    
    [Mm,Vx,Vy] = perform_nl_means(M, options);
    
    % perform partial update
    M = M*(1-lambda) + Mm*lambda;
    
    Theta = theta_list(Vx);
    Delta = delta_list(Vx);
    
    % save
    warning off;
    imwrite(clamp(Mm), [repimg name '-synth-' num2string_fixeddigit(i,2) '.png'], 'png');
    warning on;
end
fprintf('\n');


