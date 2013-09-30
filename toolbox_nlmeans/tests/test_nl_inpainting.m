% test for inpainting using NL filtering
%
%   Copyright (c) 2007 Gabriel Peyre


path(path, '../images/');
options.null = 0;

% size of input
n = 128;
n0 = [];

if ~exist('name')
    name = 'dunes';
    name = 'barb';
    name = 'tomatoes';
    name = 'lena';
    name = 'corral';
    name = 'group-people';
    name = 'reptilskin';
end

if strcmp(name, 'group-people')
    n0 = 150;
end

M = load_image(name, n0);
s = size(M,3);
M = rescale( crop(M,n) );

%% compute the mask
if not(exist('grab_mode'))
    grab_mode = 'line';
end
options.mode = grab_mode;
if not(exist('grab_radius'))
    grab_radius = 4;
end
switch grab_mode
    case 'points'
        if ~exist('U')
            options.r = grab_radius;
            U = grab_inpainting_mask(M, options);
        end
    case 'line'
        options.r = grab_radius;
        [U,options.point_list] = grab_inpainting_mask(M, options);
end
Iin = find(U(:,:,1)==Inf);
Iout = find(U(:,:,1)~=Inf);
m1 = length(Iin);


rep = 'results/inpainting/';
repimg = [rep name '/rad' num2str(grab_radius) '/'];
if ~exist(rep)
    mkdir(rep);
end
if ~exist(repimg)
    mkdir(repimg);
end
name1 = [name '-rad' num2str(grab_radius)];

B = ones(n); B(Iin) = 0;
h = ones(3); h = h/sum(h(:));
B = perform_convolution(B, h);
[D,L] = eucdist2(logical(B>=0.99));
[MapVx,MapVy] = ind2sub([n n],L);

%% initialize with random noise
M1 = M;
for i=1:s
    M1a = M1(:,:,i);
    v = randn(length(Iin),1) * std(M1a(Iout)) + mean(M1a(Iout));
    M1a(Iin) = clamp(v,min(M1a(Iout)),max(M1a(Iout)));
    M1(:,:,i) = M1a;
end

%% options of NL means
options.T = 0.06;       % width of the gaussian, relative to max(M(:))  (=1 here)
options.max_dist = 15;  % distance for search
options.ndims = 20; % number of dimension for matching (PCA) 30 is fine
options.do_patchwise = 1;

% random initial position for best match
Vx = floor( rand(n)*(n-1) ) + 1;
Vy = floor( rand(n)*(n-1) ) + 1;
% project the point outside the mask
I = sub2ind([n n], Vx,Vy);
Vx = MapVx(I); Vy = MapVy(I);
options.Vx = Vx; options.Vy = Vy;


%% options for synthesis
niter = 10;
k_schedule = [5 4 3];
options.wmax = max(k_schedule);
Tmax = 0.03; Tmin = 0.005;
use_k_schedule = 1;
use_thresh_decay = 1;
use_spacial_histo = 1;
use_exact_match = 1;

warning off;
imwrite(clamp(M), [repimg name1 '-full.png'], 'png');
imwrite(clamp(M1), [repimg name1 '-iter00.png'], 'png');
warning on;


%% iterations
kold = -1;
for i=1:niter
    progressbar(i,niter);
    
    % set width of the windows
    if use_thresh_decay
        options.T = Tmax - (i-1)/(niter-1)*(Tmax-Tmin);
    end
    % set windows half width
    ik = floor( (i-1)/niter*length(k_schedule) )+1;
    options.k = k_schedule(ik);
    if use_exact_match && i>niter*2/3
        options.T = 1e-20;
    end

    if options.k~=kold
        % compute the projection PCA matrix
        options.P = []; options.Psi = [];
        [tmp,options.P, options.Psi] = perform_lowdim_embedding(M,options);
        % input image for patch
        Ma = M1;
        for t=1:s
            Mat = Ma(:,:,t);
            Mat(Iin) = 1e6; % avoid matching
            Ma(:,:,t) = Mat;
        end
        options.Ha = perform_lowdim_embedding(Ma,options);
    end
    kold = options.k;
   
    % perform filtering
    options.Ma = M1;
    [MM,Vx,Vy] = perform_nl_means(M1, options);
    % project the point outside the boundary
    I = sub2ind([n n], Vx,Vy);
    Vx = MapVx(I); Vy = MapVy(I);
    options.Vx = Vx; options.Vy = Vy;
    if use_spacial_histo
        for t=1:s
            Mt = M(:,:,t);
            MM(:,:,t) = perform_histogram_equalization(MM(:,:,t), Mt(Iout));
        end
    end
    % impose boundary conditions
    for t=1:s
        M1a = M1(:,:,t);
        MMa = MM(:,:,t);
        M1a(Iin) = MMa(Iin);
        M1(:,:,t) = M1a;
    end

    % save current iteration
    warning off;
    imwrite(clamp(M1), [repimg name1 '-iter' num2string_fixeddigit(i,2) '.png'], 'png');
    warning on;
end
fprintf('\n');

M0 = M; 
for i=1:s
    M0a = M0(:,:,i);
    M0a(Iin) = 0;
    M0(:,:,i) = M0a;
end

%% display
clf;
ax(1) = subplot(1,2,1);
imagesc(M0); axis image; axis off;
ax(2) = subplot(1,2,2);
imagesc(M1); axis image; axis off;
colormap gray(256);
linkaxes(ax,'xy');
saveas(gcf, [rep name1 '-inpainting.png'], 'png');

warning off;
imwrite(rescale(M0), [repimg name1  '-original.png'], 'png');
imwrite(rescale(M1), [repimg name1 '-final.png'], 'png');
warning on;