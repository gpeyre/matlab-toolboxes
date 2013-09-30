ccc% Test for MCA decomposition in a wavelet+DCT basis.
% Can handle only image separation or separation+inpainting

path(path, 'toolbox/');

%% general options
do_inpainting = 1;
mask_type = 'user';
mask_type = 'chessboard';
mask_type = 'rand';

name = 'barb';
n = 128;
n0 = [];

M0 = load_image(name,n0);
M0 = rescale( crop(M0,n) ); 
M0 = rescale(sum(M0,3) );

rep = ['results/wav-dct/' name '/'];
if not(exist(rep))
    mkdir(rep);
end

%% inpainting mask computation
if do_inpainting
    if not(exist('U'))
        switch mask_type
            case 'user'
                U = grab_inpainting_mask(M0, 6);
            case 'chessboard'
                a = 5; % width of the holes
                [Y,X] = meshgrid(0:n-1,0:n-1);
                V = mod( floor(X/a)+floor(Y/a), 2);
                U = M0; U(V==0) = Inf;
            case 'rand'
                a = 60; % percent of removed pixels
                sel = randperm(n^2); 
                sel = sel(1:round(a/100*n^2));
                U = M0; U(sel) = Inf;
            otherwise
                error('Unknown mask');
        end
    end
else
    U = [];
end


% options common to all the transforms
opt.Jmin = 3;   % minimum scale for wavelets
opt.n = n;      % size of the image

components = {};
%% load the wavelet dictionary
opt.threshold_factor = 1; % reduce influency if <1
clear cpt;
cpt.options = opt;
cpt.callback =  @callback_atrou;
cpt.name = 'wav';
cpt.tv_correction = 1; % add TV penalty
components{end+1} = cpt;
%% load the local DCT dictionary
opt.w = 32;      % size of patches
opt.q = opt.w/4; % controls overlap
opt.threshold_factor = 1; % reduce influency if <1
clear cpt;
opt.dct_type = 'orthogonal4';
opt.dct_type = 'redundant';
opt.remove_lowfreq = 1; % force to 0 the lowfrequency
cpt.options = opt;
cpt.callback =  @callback_localdct;
cpt.name = 'dct';
components{end+1} = cpt;

%% options for MCA
options.niter = 60;
options.Tmax = 2.5;
options.Tmin = 0;
options.n = n;
options.InpaintingMask = U;

%% perform the mca
M = M0;
if do_inpainting
    M(U==Inf) = 0;
end
disp('--> Performing MCA.');
ML = perform_mca(M, components, options);
% recovered image
M1 = sum(ML,3);

warning off;
imwrite(rescale(M), [rep name '-original.png'], 'png');
warning on;

s = size(ML,3);
clf;
if do_inpainting
    s1 = 2;
    s2 = s;
    imageplot(M, 'Original', s1,s2,1); 
    imageplot(M1, 'Inpainted', s1,s2,2); 
else
    s1 = 1;
    s2 = s+1;
    imageplot(M1, 'Original', s1,s2,1); 
end
for i=1:s
    cpt = components{i};
    imageplot(ML(:,:,i), ['Layer ' cpt.name], s1,s2,i+s1);
    warning off;
    imwrite(rescale(ML(:,:,i)), [rep name '-' cpt.name '.png'], 'png');
    warning on;
end
saveas(gcf, [rep name '-wav-dct.png'], 'png');
