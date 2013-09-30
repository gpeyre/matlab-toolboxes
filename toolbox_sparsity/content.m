%% Toolbox Sparsity - A toolbox for sparse coding and sparse regularization
%
% Copyright (c) 2008 Gabriel Peyre
%

%% 
% The toolbox can be downloaded from Matlab Central
% http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=16204&objectType=FILE

%%
% We first includes in the path some additional useful scripts.
path(path, 'toolbox/');

%% Overview
% 
% This toolbox implements several algorithms to compute sparse expansion in redundant dictionaries and to solve inverse problems with sparse regularization 
% (and also TV regularization). 
%

%% 
% Sparse expansion of a signal |f| of size |n| in a redundant dictionary |D| (which is a matrix of size |n x p| with |p>n|, or an implicit operator)
% corresponds to the computation of coefficients |x| of size |p| such that |f=D*x| (or approximate equality |norm(y-D*x)<epsilon|). Since many |x| are possible, one 
% can assume that |x| is sparse, for instance supposing that it as a small L1 norm |sum(abs(x))|. This can be achieved by Lagrangian optimization
%

%%
% |min_x 1/2*norm(f-D*x)^2 + T*sum(abs(x))|  (1)

%%
% where |T| is increased if |epsilon| is increased. 

%%
% Sparse regularization for the resolution of |y=U*f + noise| (where |U| is a rank defficientis equivalent) corresponds to the computation of a solution 
% |f| that is sparse in a dictionary |D|, thus being |f=D*x| with small L1 norm. This is performed by solving

%%
% |min_x 1/2*norm(y-U*D*x)^2 + T*sum(abs(x))|   (2)

%%
% And the solution is computed from the optimal coefficients |x| as |f=D*x|.
%
% Problems (1) and (2) are equivalent if one replaces the dictionary |D| of (1) by the transformed dictionary |U*D| in (2). One can thus use the same 
% software to solve both (2) and (1).

%%
% TV regularization corresponds to replacing (2) by 

%%
% |min_f 1/2*norm(y-U*f)^2 + T*TV(f)|   (3)
%
% where |TV(f)| is the discrete TV norm of the signal/image |f|. Intuitively, this corresponds to assuming that the gradient of |f| is sparse, so (3) is 
% quite close to (2) if |D| is chosen a translation invariant Haar wavelet dictionary. But optimizing (3) and (2) necessitate different algorithms.

%%
% 
% * To solve (1) or (2) with |epsilon=0| (no noise), which corresponds to the minimum L1 norm solution to |f=D*x| or |y=U*D*x|, this toolbox implements the Douglas-Rachford iterative algorithm.
% * To solve (1) or (2), this toolbox implements the iterative soft thresholding algorithm. If |T| is not known but |epsilon| is known, then the implementation automatically computes the correct |T| during the iterations.
% * To solve (3), this toolbox implements Chambolle's JMIV 2002 algorithm. If |T| is not known but |epsilon| is known, then the implementation automatically computes the correct |T| during the iterations.

%% Compressed sensing of exactely sparse signals with Fourier measurements.

%% 
% First, we create an |s|-sparse signal with random spikes locations and
% random coefficients values.

% dimensionality
n = 1024;
% sparsity
s = 25; 
% generate random spikes
x0 = compute_rand_sparse(n, s, 'uniform');

%%
% We select an operator |U|, which is a pre-defined implicit operator.
% This is done by filling in a structure |options| with the correct
% parameters. Many other measurement matrices are available for Compressive
% Sensing (Gaussian, Bernouilli, Fourier, Sinus, Hadamard, etc).

% number of measurements
p = 100;
% type of matrix
clear options;
options.cs_type = 'fourier';
options.n = n;
options.p = p;

%% 
% We perform the measurements (without noise) by applying the sensing
% operator.

y = callback_sensing_rand(x0, +1, options);


%%
% We can display the signal vector and the measurements

clf;
subplot(2,1,1);
plot_sparse_diracs(x); title('Signal x_0');
subplot(2,1,2);
h = plot(real(y)); axis tight;
set(h, 'LineWidth', 2); 
set(gca, 'FontSize', 20);
title('Measurements y=U*x_0');

%%
% Since there is no noise, the recovery can be performed using
% Douglas-Rachford algorithm. Special thanks to Jalal Fadili for his help
% on this algorithm. 

% set to 1 these options if you want cool displays
options.verb = 0;
options.draw_iter = 0;
% parameter of the algorithm
options.niter = 3000;
options.x = zeros(n,1);
[xlun,lun] = perform_douglas_rachford(@callback_sensing_rand, y, options);
xlun = real(xlun);

%%
% We can display the result. This is a case of perfect recovery, |xlun=x0|!
% However, if you increase slightly the value of |s|, then the recovery
% fails ...

clf;
options.title = {'Signal' 'Recovery'};
options.val = .01; % remove small residual
plot_sparse_diracs({x0 xlun}, options);

%%
% We can also display the decrease of the L1 norm of the solution during
% the iterations, to check that the algorithm really does a good job.

clf;
h = plot(lun); axis tight;
set(h, 'LineWidth', 2);
set(gca, 'FontSize', 20);
xlabel('#iteration'); ylabel('L_1 norm |x|_1');

%% Compressed sensing of a noisy sparse signals with Gaussian measurements.

%%
% To be more robust to noise, we need to decrease the number of
% measurements. Also, we are going to use Gaussian real measurements, which
% in fact perform twice less measurements than complex Fourier
% measurements.

% sparsity
s = 20; 
% generate random spikes
x0 = compute_rand_sparse(n, s, 'uniform');

%% 
% We perform noisy measurements.
% We also select an other kind of matrix, given in explicit form (an array |U|
% of double). This is possible only for small size problems.

% Generate a matrix for sensing
U = randn(p,n);
% noiseless measurements
y = U*x0;
% make some noise
sigma = .06*std(y(:)); % noise level
y = y + sigma*randn(size(y));

%%
% Since there is some noise, we need to use iterative thresholding.
% The regularization parameter is adapted automatically to fit the noise
% level |sigma|. 

% amplification is >1 to further sparsify the solution
amplification = 2;
% norm of the noise |epsilon|
options.etgt = amplification * sigma*sqrt(p);
% initial regularization parameter
options.T = 1;
% perform iterative thresholding
options.niter = 4000;
[xlun,err,lun,lambda] = perform_iterative_thresholding( U,y,options);


%% 
% We can display the evolution of the threshold |lambda| during the iterations
% and also display the decay of the Lagrangian to see what the algorithm is
% doing.

clf;
subplot(2,1,1);
h = plot(lambda); axis tight;
set(h, 'LineWidth', 2);
set(gca, 'FontSize', 20);
title('Evolution of \lambda');
subplot(2,1,2);
lagr = .5*err.^2 + lambda(end)*lun;
h = plot(lagr); axis tight;
set(h, 'LineWidth', 2);
set(gca, 'FontSize', 20);
title('Evolution of |y-U x|^2+\lambda |x|_1');
axis([1 options.niter min(lagr) min(lagr)*1.5]);

%%
% The solution of the L1 minimization is biased, because L1 tends to
% under-estimate the value of the true coefficients (this is very much
% similar to the issue with soft thresholding). To remediate to this
% issue, it is possible to "debiased" the solution by extracting the
% support and then performing an L2 best fit (which is not biased) of the
% measurements. Be aware that this might give bad results if the
% support is baddly estimated.

% remove small amplitude coefficients because they are not very reliable
xdeb = xlun .* (abs(xlun)>.05);
% perform L2 regression
xdeb = perform_debiasing(U,xdeb,y);

%%
% We can display the result. There is no more perfect recovery because of
% the noise, but the solution is close to the original signal.


clf;
options.title = {'Signal' 'Recovery' 'Debiased'};
options.val = .01; % remove small residual
plot_sparse_diracs({x0 xlun xdeb}, options);


%% Tomography inversion of an image with TV regularization.
% We can recover from partial Fourier measurements using TV regularization.
% We add some noise to the measurement, and perform the recovery using
% Chambolle's algorithm for TV regularization.

%% 
% First we load the image.

n = 256;
M = load_image('phantom', n);
M = rescale(M,.05,.95);

%%
% Then we load the tomography mask.

nrays = 18; % Candes/Romberg/Tao set up
mask = compute_tomography_mask(n,nrays);
% number of measurements
p = sum(mask(:)==1);
% display
imageplot({M mask}, {'Signal' 'Frequencies'});

%%
% We set up |options| for the tomography operator (sub-sampling of the
% frequencies) and perform measurements.

options.size = size(M);
options.mask = mask;
y = callback_tomography(M,+1,options);
% add some noise
sigma = .01*std(y(:));
y = y + sigma*randn(size(y));

%%
% We can perform L2 pseudo-inversion (zero padding in Fourier) by just applying
% the transpose operator.

M2 = callback_tomography(y,-1,options); M2 = real(M2);
imageplot({M clamp(M2)}, {'Original' 'L2 inversion'});


%%
% We can perform inversion with TV regularization.

% set up parameter for the regularization
options.TVoperator = @callback_tomography;
options.TVrhs = y;
options.etgt = 4*sigma*sqrt(p); % target noise level
options.niter = 300;
options.niter_inner = 30; % iteration for inner loop of Chambolle's algorithm.
options.tol = 1e-4;
options.display = 0; % set to 1 for sexy display
options.verb = 0;
options.tau = 1/8;
% perform the regularization
[Mtv,err,tv,lambda] = perform_tv_denoising(zeros(n),options);

%%
% We can display the results.

clf;
imageplot({M clamp(Mtv)}, {'Original' 'Recovery'});

%%
% We can display the evolution of the regularization
% parameter and the TV energy.

clf;
subplot(2,1,1);
h = plot(lambda); axis tight;
set(h, 'LineWidth', 2);
set(gca, 'FontSize', 20);
title('Evolution of \lambda');
subplot(2,1,2);
lagr = .5*err.^2 + lambda(end)*tv;
h = plot(lagr); axis tight;
set(h, 'LineWidth', 2);
set(gca, 'FontSize', 20);
title('Evolution of |y-U x|^2+\lambda |x|_1');
axis([1 options.niter min(lagr) min(lagr)*1.5]);


%% Sparse spikes deconvolution with L1 regularization and Matching Pursuit
% sparse spikes deconvolution corresponds to inverting a band pass filter.
% This is very much used in seismic imaging where the filter is called a
% wavelet. Since the signal (reflecticity of the ground) is composed of
% isolated spikes, one uses L1 minimisation in the diract basis to recover
% the signal. This can also be performed by other sparse optimization
% technics such as greedy matching pursuit. Note that the operator |U| is
% equal to the dictionary |D| in this case, and is composed of translated
% copy of the filter. Since this results in a highly ill posed operator,
% one needs to sub-sample the dictionary a little.

%% 
% First we load a filter which is a second derivative of a Gaussian

clear options;
n = 1024;
options.sigma = .02;
h = compute_sparse_spike_filter('dergauss',n,options);
h = h/max(h);

%%
% We can display the filter and its fourier transform (to see how much
% frequencies are removed by the filtering).

hf = real(fft(h)); hf = hf / max(hf);
clf;
subplot(2,1,1);
hh = plot(fftshift(h)); axis tight;
set(hh, 'LineWidth', 2);
set(gca, 'FontSize', 20);
title('Filter');
subplot(2,1,2);
hh = plot(fftshift(hf)); axis tight;
set(hh, 'LineWidth', 2);
set(gca, 'FontSize', 20);
title('Fourier Transform');

%%
% Then we build a dictionary of sub-sampled translated copies of |h|.

q = 2; % sub-sampling (distance between wavelets)
p = n/q;
[Y,X] = meshgrid(1:q:n,1:n);
D = h(mod(X-Y,n)+1);

%%
% We can display a sub-set of elements of the dictionary.

clf;
sel = round( linspace(1,p,5) ); sel(end) = [];
hh = plot(D(:,sel)); axis tight;
set(hh, 'LineWidth', 2);
set(gca, 'FontSize', 20);
title('Translated wavelets');

%%
% Joel Tropp and Jean-Jacques Fuchs have derived in several IEEE Tr. Info.
% Theory papers some powerfull bounds of the recovery of signals with L1
% minimization. This includes the Fuchs (which depends on the sign of the coefficients), 
% ERC and WERC (which depend only on the support of the coefficients) criterions.
% We display here these three criterions to see wether a Dirac spikes train
% of spacing |Delta| between two consecutive spikes can be recovered stably
% by both matching 
% pursuit and L1 minimization (basis pursuit).
% The dashed vertical lines display the critical spacing as detected by
% these three criterions.

options.display = 1;
options.delta_max = 100;
options.verb = 0;
[minscale,crit,lgd] = compute_minimum_scale(D, options);
set(gca, 'FontSize', 20);
xlabel('\Delta');

%%
% We create a sparse spikes train, with a default in the middle (two close
% spikes). Note that we select the spacing between the spikes slightly below the minimum scale detected
% by ERC<1, which causes sparsity methods to become unstable, and matching pursuit to fail.

options.delta = 23; % regular spacing
options.delta0 = 3; % create default that is not recovered by BP
x = compute_rand_sparse(p, 0, 'seismic', options);

%%
% The signal is the convolution with the filters plus noise.

sigma = .03;
y = D*x + randn(n,1)*sigma;
% display
clf;
subplot(2,1,1);
plot_sparse_diracs(x);
title('Sparse spikes x');
subplot(2,1,2);
hh = plot(y); axis tight;
set(hh, 'LineWidth', 2);
set(gca, 'FontSize', 20);
title('Signal y=Dx+w');

%%
% We can recover the signal using L1 minimisation plus debiasing.

options.T = 1;
options.niter = 1000;
options.etgt = 5*sigma*sqrt(n); % noise level
options.x = zeros(p,1); % initial solution
[xlun,err,lun,Tlist]  = perform_iterative_thresholding(D,y,options);
% perform debiasins
xlun(abs(xlun)<.08)=0;
xdeb = perform_debiasing(D,xlun,y);
% display
options.title = {'Signal' 'L1' 'Debiased'};
plot_sparse_diracs({x xlun xdeb}, options);


%%
% Another way to solve the problem is using greedy matching pursuit
% procedure. Most of the time, it gives results comparable or a little worse
% than basis pursuit, but it is faster. Here, we are in a critical set up
% where basis pursuit performs better than matching pursuit. The toolbox
% implements basis matching pursuit and its orthogonalized version. Only
% the basic pursuit is implemented with implicit operator, and orthogonalized
% pursuit requires the computation of many pseudo inverse (so it is quite slow).

% parameters for the pursuits
options.tol = 3*norm(w)/norm(y);
options.nbr_max_atoms = 80;
xmp = perform_mp(D,y,options);
xomp = perform_omp(D,y,options);
% display
options.title = {'Signal' 'Matching pursuit' 'Orthogonal Matching Pursuit'};
plot_sparse_diracs({x xmp xomp}, options);


%% Image inpainting with sparse wavelets regularization
% Image inpainting consists in filling missing pixels in an image. Many inpainting
% algorithm rely on a diffusion PDE that try to propagate information from
% the boundary of the missing region. An alternative methods consist in
% treating this problem as an inverse problem. Inpainting indeed
% corresponds to a diagonal operator |U| with 0 and 1 on the diagonal (0
% indicating missing pixels). This can be regularized using iterative
% thresholding. Since there is no noise, we want to use a very low
% threshold. This can be achieved by decaying the threshold during the
% iterations, and corresponds to the Morphological Component Analysis
% strategy of Jean-Luc Starck, Miki Elad and David Donoho.
% See also the homepage of Jalal Fadili for many example of inpainting.

%%
% First we load an image and a mask.

clear options;
n = 128;
M = load_image('lena');
M = rescale( crop(M,n) );
options.rho = .7; % remove 70% of pixels
mask = compute_inpainting_mask('rand',n,options);
options.mask = mask;

%%
% Apply the operator to remove pixels.

y = callback_inpainting(M, +1, options);
imageplot({M y}, {'Image' 'Data to inpaint'});


%% 
% First we set up the dictionary: a translation invariant wavelet
% dictionary.

options.Jmin = 4;
options.wavelet_type = 'biorthogonal';
options.wavelet_vm = 3;
options.D = @callback_atrou;

%% 
% Then we set up decaying thresholding parameter for the MCA solver.

options.thresh = 'soft';
options.Tmax = .1;
options.Tmin = 0;
options.niter = 200;
options.tau = 1;
options.drawiter = 0;
options.verb = 0;

%%
% Perform the inpainting using sparse wavelet expansion.

% do the resolution
[MW,err,lun,Tlist] = perform_iterative_thresholding(@callback_inpainting, y, options);
% retrieve the image from its coefficients
Mlun = callback_atrou(MW, +1, options);

%% 
% Display the result

imageplot({M clamp(Mlun)},{'Original' 'Recovered'});

%% Dictionary Learning
% Instead of using a fixed dictionary |D| to compress or denoise a signal
% with sparse coding (or to solve an inverse problem), it is possible to
% optimize |D| in order to sparsify a set of given exemplar. In practice,
% these exemplar are taken as small patches extracted from tons of natural
% images. As first shown by Olshausen and Field, the optimized atoms then
% ressemble those of an oriented (steerable) wavelet frame.

%%
% The parameter of the dictionary

% size of the patches
w = 12; n = w^2;
% redundancy of the dictionary
redun = 1.5;
% number of atoms
p = round( n*redun );
% overtraining factor
overtraining = 6;
% number of example
m = round( overtraining*p );

%% 
% First we load a set of images
Mlist = load_image({'lena' 'barb' 'boat'});

%%
% Then we extract a large set of patches from these images

Y = load_random_patches(Mlist,w,m);


%% 
% Parameters for the learning of the dictionary (optimization of the atoms
% using an iterative method).

% number of atoms
options.K = p;
% number of iterations
options.niter_learning = 60;
% solver used for the sparse coding stage
options.sparse_coder = 'omp';
options.sparse_coder = 'mp';
% sparsity targeter for the sparse code
options.nbr_max_atoms = 5;
% initialization method
options.options.init_dico = 'input';

%%
% Perform the learning using an optimization procedure called the MOD
% (Method of Direction) algorithm.

options.learning_method = 'mod';
options.verb = 0;
[D1,X1,E1] = perform_dictionary_learning(Y,options);

%%
% Display the vectors of the dictionary as small images.
% You can notice that the texture of barbara is occupating many atoms !

clf;
options.ndim = 2;
options.normalization = 'clamp';
display_dictionnary(D1, X1, [10 14], options );