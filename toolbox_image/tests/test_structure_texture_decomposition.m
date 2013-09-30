% test for TV-L2 (Osher-Rudin-Fatemi) decomposition
% and TV-Gabor (Aujol-Gilboa) decomposition

n = 128;

path(path, 'toolbox/');

name = 'barb';

% to save results
rep = ['results/separation/' name '/'];
if not(exist(rep))
    mkdir(rep);
end

if strcmp(name, 'barb')
    n0 = [];
else
    n0 = n;
end

M = load_image(name, n0);
M = rescale( sum( crop(M,n),3) );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TV-L2 decomposition
options.dt = 0.1249;
lam = 1/20 * 215/max(M(:));
options.niter = 300;
options.use_gabor = -1; % 1:TV-Hilbert, -1:TV-L2
lam = 1/50 * 215/max(M(:));

options.p0x = []; options.p0y = [];
disp('--> Performing TV-L2 decomposition.');
options.p0x = []; options.p0y = [];
[u0,v0,options.p0x,options.p0y] = perform_tv_hilbert_projection(M,[],lam,options);

clf;
imageplot(u0,'cartoon part',1,2,1);
imageplot(v0,'textured part', 1,2,2);

warning off;
imwrite( clamp(M), [rep name '-original.png'] );
imwrite( clamp(u0), [rep name '-tvl2-structure.png'] );
imwrite( rescale(v0), [rep name '-tvl2-texture.png'] );
warning off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TV-Gabor decomposition
lam = 0.005 * 215/max(M(:)); %parameter for TV-Hilbert
options.niter = 1000;  % number of iterations for TV-Hilbert
options.dt=0.001; %step time for TV-Hilbert
options.use_gabor = 1; % 1:TV-Hilbert, -1:TV-L2

disp('--> Performing TV-Gabor decomposition.');
[u0,v0,options.p0x,options.p0y] = perform_tv_hilbert_projection(M,@compute_local_fourier_operator,lam,options);

clf;
imageplot(u0,'cartoon part',1,2,1);
imageplot(v0,'textured part', 1,2,2);

warning off;
imwrite( clamp(u0), [rep name '-tvk2-structure.png'] );
imwrite( rescale(v0), [rep name '-tvk2-texture.png'] );
warning off;