% test for tensor voting
% Use the tensor structure direction field as an input.

sigma = 0.5;
mu = 3;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a test image
disp('--> Building test image.');
Radius = 0.9;

n = 75;
x = -1:2/(n-1):1;
y = linspace(0,1,n);
[X,Y] = meshgrid(x,x);

M = double( abs( X.^2+Y.^2-Radius^2 )<0.03 );
M = double( X.^2+Y.^2<Radius^2 );

% add some noise in the middle of the image
sigma = 0.8;    % level of the noise
r = 0.05;        % width of the gap
I = find(abs(X-0.5)<r);
M = M + sigma*randn(n).*(abs(X-0.5)<r);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% structure tensor computation
disp('--> Compute tensor structure.');
h = compute_gaussian_filter([11 11],0.6/n,[n n]);
Mh = perform_convolution(M,h);
ST = compute_rigidity_tensor(Mh);
h = compute_gaussian_filter([11 11],0.3/n,[n n]);
for i=1:4   % convole each entry
    ST(:,:,i) = perform_convolution( ST(:,:,i) ,h);
end
[e1,e2,l1,l2] = perform_tensor_decomp(ST);   % e1 is a regularized gradient

% As the main direction is the direction orthogonal to the discontinuity 
% we must swap the eigenvectors.
T = perform_tensor_recomp(e2,e1,l1,l2);

clf;
subplot(2,2,1);
plot_vf(e2,M);
axis off;
colormap gray(256);
title('Original Image and Field');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% an example of field
A = compute_stick_tf([1,0.3],n);
[e1,e2,l1,l2] = perform_tensor_decomp(A); 

subplot(2,2,2);
plot_vf(e1,l1);
axis off;
title('Example of Stick Field');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% peform voting
disp('--> Performing tensor voting.');
sigma = 0.2;
c = 0.1*sigma;
T = perform_voting(T,sigma,c);
disp('--> Performing tensor decomposition.');
[e1,e2,l1,l2] = perform_tensor_decomp(T,'abs');

% re-orient vf
e1 = perform_vf_reorientation(e1);
e2 = perform_vf_reorientation(e2);

% plot e1
subplot(2,2,3);
plot_vf(e1,l1);
axis off;
title('New Field with Eigenvalue');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% curve extraction via zeros crossing of A=ds/dn where s=l1-l2, n=e2,
% orthogonal to e1
s = l1-l2;
g = compute_grad(s);
A = prod_vf_vf(g,e2);

% extract only the 20 longest curves
options.min_length = 20;
c_list = perform_curve_extraction( A,0,10 );


% plot curves
subplot(2,2,4);
plot_curve(c_list,A);
axis off;
title('Extracted Curves');
saveas(gcf, 'test_structure', 'png');