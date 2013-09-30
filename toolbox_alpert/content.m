%% Toolbox Alpert Transform - A toolbox to perform the Fast Multiwavelet Alpert Transform
%
% Copyright (c) 2008 Gabriel Peyre
%

%% 
% The toolbox can be downloaded from Matlab Central
% http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=5643&objectType=FILE

%% Setting-up the path and compiling the mex files.

%%
% We add some directories to the path.
path(path, 'toolbox/');


%% Introduction to the Alpert Multiwavelets
% The Alpert transform is a multiwavelets transform based on orthogonal polynomials. 
% It was originally designed for the resolution of partial differential 
% and integral equations, since it avoids boundary artifact and can be 
% used with an arbitrary sampling.
 
%%
% The reference for the numerical algorithm:

%%
% * Bradley K. Alpert, _Wavelets and Other Bases for Fast Numerical Linear Algebra_, in Wavelets: A Tutorial in Theory and Applications,  C. K. Chui, editor, Academic Press, New York, 1992.

%%
% And more theoretical (continuous setting):

%%
% * B.K. Alpert, _A class of bases in L^2 for the sparse representatiion of integral operators_, in SIAM J. Math. Anal., 24 (1993), 246-262.
 
%%
% The strengh of this transform is that you can transform data sampled irregularly. 
% Of course this algorithm runs in linear time, i.e. O(n). The use of multiwavelets 
% remove any boundary artifact (which are common with wavelet of support > 1, e.g.
% Daubechies wavelets), but the price to pay is that the wavelets functions are *not* 
% continuous, they look like the Haar basis functions. So do not use this transform to 
% compress data that will be seen by human eyes (although the reconstruction error 
% can be very low, the reconstructed function can have some ugly steps-like artifacts).
 
%%
% In this toolbox, you can transform a signal (1D, 2D, nD) with arbitrary length and 
% arbitrary sampling (you *must* each time provide sampling locations 
% in a parameter |pos|).

%%
% _The number of vanishing moments_ (which is also the degree of the polynomial 
% approximation+1) is set via the parameter 'k' (default=3). You can provide 
% different numbers of vanishing moments for X and Y axis, using |k=[kx,ky]| 
% (default |k=[3,3]|).

%%
% _Dichotomic Grouping Functions_: All the transforms require to recursively split the data. For each transform
% a default splitting rule is used, but you can provide your own 
% partition via a parameter |part|. The way you split the data
% (isotropicaly, using a prefered direction, etc) will lead to various
% behaviour for the transform (of course this is relevant for 2D transform only, 
% since the splitting rule in 1D is always the same).
% The function to perform an automatic grouping is 'dichotomic_grouping'
% (it can provide orthogonal split along the axis or isotropic using k-means).

%% 1D Alper transform at regularly sampled locations
% The Alpert transform act just like a traditional wavelet transform (see my Wavelets Toolbox for more details). It
% is an orthogonal transform, which wavelet of minimal support (just like
% the haar transform) and with as many vanishing moments as desired. The transform algorithm is O(n), which is 
% similar to the traditional wavelet transform. The
% main drawback is that the basis wavelet functions are not regular, which
% might lead to visually deapleasant artifacts after approximation or
% compression of the signal.

%%
% We first load a signal
n = 512;
f = rescale( load_signal('Piece-Regular',n) );

%%
% Then we perform the forward transform.

% sampling location
x = (0:n-1)/n;
% number of vanishing moments
k = 4;
[fw,info] = perform_alpert_transform_1d(f, x, k, +1);

%%
% To perform approximation, we discard the smallest amplitude coefficients,
% and then reconstruct.

fwT = keep_biggest(fw, round(.1*n) );
f1 = perform_alpert_transform_1d(fwT, x, k, -1);

%% 
% For the display, we reverse the coefficients using |fw(end:-1:1)|, so
% that coarse scale coefficients are located on the left. Each vertical red dashed 
% line is the separation between two wavelet scales.

clf;
subplot(3,1,1);
plot(f); axis tight; title('Signal');
subplot(3,1,2);
plot_wavelet(fw(end:-1:1), 1); axis([1 n -.3 .3]); title('Coefficients');
subplot(3,1,3);
plot(f1); axis tight; title('Approximated');

%%
% You can display the Alpert basis vectors by taking the inverse transform
% of a dirac. For each scale and each positon, there are |k| Alpert wavelets,
% this is why this basis is called a multiwavelet basis.

jlist = [2 3 4];
klist = 1:k;
clf;
for a=1:length(jlist)
    f = zeros(n,k);
    lgd = {};
    for b=1:length(klist)
        % construct a diract at the correct location
        t = klist(b); j = jlist(a);
        lgd{end+1} = ['k=' num2str(t)];
        fw = zeros(n,1);
        I = find( info.l==j & info.k==t );
        fw(I(round(length(I)/2))) = 1;
        % inverse transform
        f(:,b) = perform_alpert_transform_1d(fw, x, k, -1);
    end 
    % display
    subplot( length(jlist), 1, a);
    plot(f); axis tight;
    title(['j=' num2str(j)]); legend(lgd)
end



%% 1D Alpert transform on irrregular samples 
% The Alpert transform can be used to process signal irregularely sampled.
% The vectors keep their properties (minimal support, vanishing moments,
% orthogonality) even on this setting, which is a major advantage of this
% transform, and makes it unique among wavelet transform. 


% generate regular and random sampling
xreg = (0:n-1)/n;
xireg = rand(n,1).^3; xireg = sort(rescale(xireg));
% generate a piecewise smooth signal for this sampling
f = mod(xireg.^2,.2);
k = 3; % number of VM
fw0 = perform_alpert_transform_1d(f, xreg, k, +1);
fw1 = perform_alpert_transform_1d(f, xireg, k, +1);

%%
% A comparison of the coefficients with correct and wrong sampling reveals
% that the coefficients with correct sampling are zero away from
% singularities. This is because |f| is a piecewise second degree
% polynomial, and the Alpert basis has 3 vanishing moments.

clf; eta = .15;
subplot(3,1,1);
plot(xireg, f, '.-'); axis tight; title('Signal and sampling');
subplot(3,1,2);
plot_wavelet(fw0(end:-1:1),1); axis([1 n -eta eta]); title('Coefficient with (wrong) uniform sampling');
subplot(3,1,3);
plot_wavelet(fw1(end:-1:1),1); axis([1 n -eta eta]); title('Coefficient with (correct) non-uniform sampling');


%% 2D Alpert transform 
% Use function |perform_alpert_transform_2d|.
% The default spliting rule for this transform use
% the X axis (but you can change). So by default, it is not a fully 2D wavelet 
% construction (the transform is pyramidal only on 1 dimension, on  the other 
% this is just a fixed polynomial basis  with a given number of slices).
% But if you provide your own subdivision (via parameters |part| for user-defined
% or |part_type| for automatic)
% then the transform can become isotropic (for example
% set |part_type='kmeans'|).


%%
% Generate a random 2D sampling
n = 20000;
pos = rand(2,n);

%%
% Compute a 2D signal by evaluating Lena image at the sampling locations
n0 = 256;
M = load_image('lena');
M = rescale(crop(M,n0));
x = linspace(0,1,n0);
[X,Y] = meshgrid(x,x);
f = interp2(X,Y,M,pos(2,:),pos(1,:));

%%
% Display the image, its non-uniform sampling, and the linear
% interpolation.

% linear interpolation
F  = griddata(pos(2,:),pos(1,:),f, X,Y);
% display
clf;
subplot(1,2,1);
hold on;
imagesc(x,x,M); axis image; axis off; axis ij;
title('Image and sampling (sub-set)');
plot(pos(1,1:1000),pos(2,1:1000), 'r.'); 
hold off;
subplot(1,2,2);
imageplot(F, 'Interpolated');

%%
% Compute the 2D Alpert transform.

k = [2 2]; % vanishing moments
options.part_type = '2axis'; % partitionning type
[fw,info] = perform_alpert_transform_2d(f, pos, k, +1, options);


%%
% Approximate the coefficients by thresholding, and display on a regular
% grid the data using interpolation.

% approximate
fwT = keep_biggest(fw, round(.1*n) );
% inverse transform
[f1,info] = perform_alpert_transform_2d(fwT, pos, k, -1, options);
% interpolate
F1 = griddata(pos(2,:),pos(1,:),f1,X,Y);
% display
clf;
imageplot({F clamp(F1)}, {'Original' 'Approximated'});


%% nD and more exotic transforms
% Use the function 'perform_alpert_transform_nd' to transform
% data in arbitrary dimension.
% The clustering uses function 'dichotomic_grouping'.
% This transform can also be used to transform 
% data lying on a manifold embedded in R^d.
% For instance, the function 'test_spherical'
% gives an example which transform data lying on a sphere.

%%
% We first make a random sampling on the sphere and compute a function
% defined on this sampling.

n = 2048;
k = 3;
pos = randn(3,n);
pos = pos ./ repmat( sqrt(sum(pos.^2,1)), [3 1]);
f = rescale( pos(1,:).^2 - pos(2,:).*pos(3,:) );
f = cos(10*f(:));
% add a discontinuity to make the approximation problem harder
f = rescale(f - rescale(pos(3,:)'>0) );

%%
% Then we perform iterative clustering to define hierachical sets of
% embedded points on the sphere.


clear options;
options.ptmax = k^2;
options.part_type = 'kmeans';
[part,B,G] = dichotomic_grouping(pos,options);
% display
clf;
subplot(1,2,1);
plot_spherical_partition(pos,part);
title('Clustering of sampling locations');
subplot(1,2,2);
plot_spherical_function(pos,f);
title('Function to transform (interpolated)');

%%
% Now we can perform the Alpert transform of the function using this
% specific spherical sampling.

[fw,info,part] = perform_alpert_transform_nd(f, pos, k, 1, options);

%%
% Approximation is perform by simply thresholding the set of coefficients
% to keep only largest amplitude coefficients.

fwT = keep_biggest(fw, round(n*.08));
options.part = part;
f1 = perform_alpert_transform_nd(fwT,pos,k, -1,options);

%%
% At last we can display the approximated function. Since the original
% function is smooth, the approximation produces a small error.

clf;
subplot(1,2,1);
plot_spherical_function(pos,f);
title('Original signal');
subplot(1,2,2);
plot_spherical_function(pos,f1);
title('Approximated signal');



%% Other features of the toolbox
% The toolbox implements several other features, among which:

%%
%
% * Functions |build_alpert_matrix| and |build_alpert_matrix_2d| to build the transformation matrix, but these are quite slow.
% * Sliced 2D transforms

%%
% A sliced tranform is implemented in |function perform_alpert_transform_2d_sliced|. 
% It divides the points into |s| slices and perform a 2D transform on 
% each slice. You have to provide a number of slices, which  fixes the 
% resolution on the Y axis (as if you were decomposing the function on 
% a wavelet basis for the X direction, and on a box spline basis for the 
% Y direction). This is not very usefull unless you want to 
% set a given resolution on the X axis.