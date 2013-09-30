% test for image warping
M = load_image('barb');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% extract eyes %%%
clf;
imagesc([0,1],[0,1],M);
title('Click on the each eye');
colormap gray(256);
axis image; axis off;
[x,y,b] = ginput(2);

% position of the eye
u = [x(1) y(1)];
v = [x(2) y(2)];

u1 = [0.3,0.3];
v1 = [0.7,0.3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% computations %%%
M1 = perform_image_similitude(M,u,u1,v,v1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% display %%%
clf;
subplot(1,2,1);
hold on;
imagesc([0,1],[0,1],M);
plot([u(1) v(1)], [u(2) v(2)], '*');
hold off;
title('Original');
colormap gray(256);
axis image; axis off; axis ij

subplot(1,2,2);
hold on;
imagesc([0,1],[0,1],M1);
plot([u1(1) v1(1)], [u1(2) v1(2)], '*');
hold on;
title('Warped');
colormap gray(256);
axis image; axis off; axis ij