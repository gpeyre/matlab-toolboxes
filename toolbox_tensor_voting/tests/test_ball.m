% test for tensor voting using only the ball tensor field.
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% a test image
disp('--> Building test image.');
Radius = 0.9;

n = 75;
x = linspace(0,1,n);
[X,Y] = meshgrid(x,x);

M = double( abs( X.^2+Y.^2-Radius^2 )<0.03 );
% M = X.^2+Y.^2<Radius^2;

% add some holes
nb_holes = 200;
I = find(M==1);
J = floor(rand(nb_holes,1)*length(I))+1;
M(I(J)) = 0;


clf;
subplot(2,2,1);
imagesc(M);
axis xy;
axis off; axis square;
colormap gray(256);
title('Original Image');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% example of field
A = compute_ball_tf(75);
[e1,e2,l1,l2] = perform_tensor_decomp(A); 
subplot(2,2,2);
plot_vf(e1,l1);
axis off;
title('Example of Ball Field');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% peform voting
disp('--> Performing tensor voting.');
sigma = 0.2;
c = 0.1*sigma;
T = perform_voting_ball(M,sigma,c);
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
% curve extraction via zeros crossing of A=ds/dn where s=l1-l2, n=e2,
% orthogonal to e1
s = l1-l2;
g = compute_grad(s);
A = prod_vf_vf(g,e2);

% extract ten curves
options.min_length = 20;
c_list = perform_curve_extraction( A,0,10 );


% plot curves
subplot(2,2,4);
plot_curve(c_list,A);
axis off;
title('Extracted Curves');


saveas(gcf, 'test_ball', 'png');