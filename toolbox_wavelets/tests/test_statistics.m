% some statistical tests


path(path, '../../images/');

if not(exist('name'))
name = 'line_vertical_blurred';
name = 'sawtooth';
name = 'line_vertical';
name = 'flow2';
name = 'flow2';
name = 'noise';
name = 'turbulence';
name = 'corral';
name = 'lena';
name = 'barb';
name = 'reptilskin';
name = 'hair';
name = 'disk';
name = 'square';
end
if not(exist('wavtype'))
wavtype = 'redun';
wavtype = 'ortho';
end

rep = ['results/stats/'];

save_images = 1;

if exist(rep)~=7
    mkdir(rep);
end

n = 512;
n0 = [];
if strcmp(name, 'disk') || strcmp(name, 'square')
    n0 = n;
end

M = load_image(name, n0); M = sum(M,3);
M = rescale(crop(M,n));
n = size(M,1);

options.wavelet_type = 'biorthogonal';
options.wavelet_vm = 4;
Jmin = 4;
Jmax = log2(n)-1;

switch wavtype
    case 'ortho'
        MW = perform_wavelet_transform(M, Jmin, +1, options);
        %% scale selection
        j = Jmax;
        [selx,sely] = compute_quadrant_selection(j,2);
        MW1 = MW(selx,sely);
        % father
        [selx,sely] = compute_quadrant_selection(j-1,2);
        MWf = MW(selx,sely);
        sel = floor( 1:0.5:nj/2+0.5 );
        MWf = MWf(sel,sel);
        % brother
        [selx,sely] = compute_quadrant_selection(j,1);
        MWb = MW(selx,sely);
    case 'redun'
        MW = perform_atrou_transform(M,Jmin,options);
        MW1 = MW{1};
        MWf = MW{4};
        MWb = MW{2};
end
nj = size(MW1,1);

if save_images
    warning off;
%    imwrite(rescale(M),[rep name '-original.png'], 'png');
    warning off;
end

%% compute 1D histogram
nb_bins = 501;
[N,X] = histo(MW1(:),nb_bins);
N = N/prod( size(MW1) );
clf;
plot( X, log2(N) );
t = max(abs(MW1(:)))*0.8;
axis( [-t t -16 max(log2(N))] );

if save_images
%    saveas(gcf, [rep name '-histo.eps'], 'epsc');
%    saveas(gcf, [rep name '-histo.png'], 'png');
end

% parameters for display
pdisp = 800;
ms = 15;
nbins = 50;

sibling_list = { 'same','same','same','same','same','father','brother' };
delta_list = { [1 0] [0 1] [2 0] [0 2] [4 0] [0 0] [0 0] };
ntests = length(sibling_list);
nrows = 3;
ncols = ceil( (ntests+2)/nrows );

clf;
subplot(nrows, ncols, 1);
imageplot( M, 'Original' );
subplot(nrows, ncols, 2);
imageplot( MW1, 'Wavelets coefs' );

for itest = 1:ntests
    sibling = sibling_list{itest};
    d = delta_list{itest};
    
    switch sibling
        case 'same'
            MW2 = MW1;
        case 'father'
            MW2 = MWf;
        case 'brother'
            MW2 = MWb;
    end
    
    %% select neighboring coefficients
    [Y,X] = meshgrid(1:nj,1:nj);
    I = find( X+d(1)>0 & X+d(1)<=nj & Y+d(2)>0 & Y+d(2)<=nj );
    v1 = MW1(I);
    J = X(I)+d(1) + (Y(I)+d(2)-1)*nj;
    v2 = MW2(J);

    %* compute covariance matrix
    w1 = v1-mean(v1);
    w2 = v2-mean(v2);
    kappa1 = 2*mean(w1.*w2)/(mean(w1.^2)+mean(w2.^2));
    w1 = w1/sqrt(sum(w1.^2));
    w2 = w2/sqrt(sum(w2.^2));
    kappa = sum(w1.*w2);
 
    % bounds
    t1 = max(abs(v1))*0.7;
    t2 = max(abs(v2))*0.7;
    t1 = std(v1)*3;
    t2 = std(v2)*3;

    
    %% display normalized pairwise
    x1 = linspace(-t1,t1,nbins)';
    x2 = linspace(-t2,t2,nbins)';
    h = compute_pairwise_histogram( cat(2,v1,v2), x1, x2, options);
    H = h./repmat( max(h,[], 1), [nbins 1] );
    H(isnan(H)) = min(H(:));
    
    str = sibling;
    if d(1)~=0 || d(2)~=0
        str = [str ',\delta=(' num2str(d(1)) ',' num2str(d(2)) ')'];
    end
    
    subplot(nrows, ncols, itest+2);
    hold on;
    imageplot(-H);
    sel = randperm(length(v1)); sel = sel(1:pdisp);
    % aa = plot( v1(sel),v2(sel), '.' );  axis([-t1 t1 -t2 t2]);
    % set(aa, 'MarkerSize', ms);    
    title([str ',c=' num2str(kappa*100, 2) '%']);
    hold off;

end
saveas(gcf, [rep name '-' wavtype '-histo-pairw.png'], 'png');