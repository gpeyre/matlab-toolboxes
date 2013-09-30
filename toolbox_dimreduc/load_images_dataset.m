function M = load_images_dataset(name, options)

% load_images_dataset - load a dataset of image
% 
%	M = load_images_dataset(name, options);
%
%	name should be either 'binaryalphadigs', 'digits', 
%       'frey_rawface' : From Brendan Frey. Almost 2000 images of Brendan's
%               face, taken from sequential frames of a small video. Size: 20x28.
%		'olivettifaces' : Grayscale faces 8 bit [0-255], a few images of several different people.
%               400 total images, 64x64 size. From the Oivetti database at ATT
%       'umist_cropped' : Grayscale faces 8 bit [0-255], a few images (views) of 20 different people.
%               575 total images, 112x92 size, manually cropped by Daniel
%               Graham at UMist.
%               Set options.nclass=1:20 to load all the data.
%       'yaleface' : The Yale Face Database (size 6.4MB) contains 165
%               grayscale images in GIF format of 15 individuals. There are 11 images per subject, 
%               one per different facial expression or configuration: center-light, w/glasses, 
%               happy, left-light, w/no glasses, normal, right-light, sad, sleepy, surprised, and wink. 
%               See : http://cvc.yale.edu/projects/yalefaces/yalefaces.html
%       'yaleface_b' :
%               See : http://cvc.yale.edu/projects/yalefacesB/yalefacesB.html
%       'teapots' : see
%               http://www.seas.upenn.edu/~kilianw/Downloads/Data.html
%       'teapots-bw' : same as teapot but in bw (simpler)
%       'edges', 'disks' : synthetic datasets.
%		(more to come)
%	options.nclass is the number of the sub-class (e.g. the number of the digit
%		for 'binaryalphadigs' database)
%	options.dim is the dimension of the images to be loaded (if not set automatically).
%	options.smoothing is to add a smoothing (width of the kernel, eg. 0.05).
%
%   For face datasets, check http://www.face-rec.org/databases/
%
%   Copyright (c) 2005 Gabriel Peyré

% name can be 'binaryalphadigs'

options.null = 0;

rep = 'data/';

if isfield(options, 'nclass')
    nclass = options.nclass;
else
    nclass = 1;
end


if isfield(options, 'nbr')
    nbr = options.nbr;
else
    nbr = 200;
end


if isfield(options, 'dim')
    dim = options.dim;
else
    dim = [32 32];
end


if isfield(options, 'smoothing')
    smoothing = options.smoothing;
else
    smoothing = 0;
end


if nargin<2
    nclass = 1;
end

switch lower(name)

    case 'binaryalphadigs'
        % this is data set of 39 binary images
        load([rep 'binaryalphadigs']);
        p = numclass;
        numclass = min(p,nclass);
        n = size(dat{1});
        M = zeros(n(1),n(2),p);
        for i=1:p
            M(:,:,i) = dat{numclass + (i-1)*36};
        end
        M = 1-M;
        
	case 'teapots'
        % this is data set of 400 color images
        M = zeros(76,101,3,400);
        for i=1:400
            M(:,:,:,i) = double( imread([rep 'teapots/teapot-' num2str(i) '.jpg']) );
        end
        M = rescale(M);
        % sub-sample if necessary
        if nbr<400
            sel = randperm(400); sel = sel(1:nbr);
            M = M(:,:,:,sel);
        end
        
    case 'teapots-bw'
        M = load_images_dataset('teapots', options);
        M = sum(M,3)/3; % bw images
        M = reshape(M, [size(M,1),size(M,2),size(M,4)]);
        m = min(size(M,1),size(M,2));
        % square images
        M = M( round(end/2-m/2+1:end/2+m/2) , round(end/2-m/2+1:end/2+m/2) ,:);
        
    case 'digits'
        rep = [rep 'digits/'];
        for i=1:nbr
            [m,map] = imread([rep 'digit' num2str(i) '.gif'], 'gif');
            map = map(:,1);
            M(:,:,i) = map(m+1);
        end

    case 'yaleface'
        nbr = min(nbr,166);
        rep = [rep 'yalefaces/'];
        n = 64;
        M = zeros(n,n,nbr);
        for i=1:nbr
            m = imread([rep 'yaleface_' num2str(i) '.png'], 'png');
            M(:,:,i) = sum(m,3);
        end

    case 'frey_rawface'
        load([rep 'frey_rawface']);
        nbr = min(nbr,size(ff,2));
        for i=1:nbr
            m = ff(:,i);
            M(:,:,i) = reshape(m,20,28);
        end
        
    case 'olivettifaces'
        load([rep 'olivettifaces']);
        nbr = min(nbr,size(faces,2));
        for i=1:nbr
            m = faces(:,i);
            M(:,:,i) = reshape(m,64,64);
        end
        
    case 'umist_cropped'
        load([rep 'umist_cropped']);
        M = [];
        for i=nclass
            M = cat(3, M, facedat{i});
        end  
        
    case 'disks'
        x = linspace( 0,1, dim(1) );
        [Y,X] = meshgrid(x,x);
        if isfield(options, 'radius')
            r = options.radius;
        else
            r = 0.15;
        end
        M = zeros(dim(1),dim(1),nbr);
        for i=1:nbr
            c = r + rand(2,1)*(1-2*r);
            M(:,:,i) = (X-c(1)).^2 + (Y-c(2)).^2 < r^2;
        end
        
        
    case 'edges'
        x = linspace( 0,1, dim(1) );
        [Y,X] = meshgrid(x,x);
        if isfield(options, 'radius')
            r = options.radius;
        else
            r = 0.15;
        end
        M = zeros(dim(1),dim(1),nbr);
        for i=1:nbr
            m = [];
            % force edge to reach 5% of total area of image
            while sum(m(:))<0.05*dim(1)^2
                t = rand(1)*2*pi;
                X = 2*X-1; Y = 2*Y-1;
                X = X + 2*rand(1)-1;
                Y = Y + 2*rand(1)-1;
                y = -sin(t)*X(:) + cos(t)*Y(:);
                m = reshape( y>0, dim(1),dim(1) );
            end
            M(:,:,i) = m;
        end
        
    otherwise
        
        error('Unkwnown dataset');


end

M = double(M);

if smoothing>0
    n = size(M,1);
    h = compute_gaussian_filter(ceil([n n]/4)*2+1,smoothing,n);
    for i=1:size(M,3)
        M(:,:,i) = conv2(M(:,:,i), h, 'same');
    end
end

