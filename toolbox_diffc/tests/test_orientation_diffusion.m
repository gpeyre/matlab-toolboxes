% test for orientation diffusion


Dt = 1.5;

m = 2*pi;

test_fs = 0;
test_prog = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show in full screen
if test_fs
n = 100;
M = rand(n)*m;

nsteps = 100;
sigma = 0.1; % in absolute coords
t = n*sigma;

M = perform_orientation_diffusion(M,t,40,m);
vf = zeros(n,n,2);
vf(:,:,1) = cos(M);
vf(:,:,2) = sin(M);
options.display_streamlines = 1;
clf;
plot_vf(vf,cos(M),options);
axis tight; axis square;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show progression
if test_prog
n = 50;
nb_iter = 6;
nbsteps = 20;
M = rand(n)*m;
sigma = 0.02; % in absolute coords
Dt = n*sigma;
options.display_streamlines = 1;

clf;
for i=1:nb_iter
    % perform iteration
    M = perform_orientation_diffusion(M,Dt*i,nbsteps,m);
    vf = zeros(n,n,2);
    vf(:,:,1) = cos(M);
    vf(:,:,2) = sin(M);
    % compute energy
    dt = 0.2;
    h = [ 0 dt/4 0; dt/4 1-dt dt/4; 0 dt/4 0 ];
    E = perform_convolution( exp( 1i * (2*pi/m) * M), h );
    E = abs(E).^2;
    % plot everything 
    subplot(2,3,min(i,6));
    plot_vf(vf,E,options);
    title(sprintf('Time %.2f',Dt*i));
    axis tight; axis square;
    colormap gray(256);
end
end