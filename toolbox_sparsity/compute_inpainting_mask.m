function mask = compute_inpainting_mask(type,n,options)

% compute_inpainting_mask - compute an inpainting mask.
%
%   mask = compute_inpainting_mask(type,n,options);
%
%   mask=1 : data removed
%
%   Copyright (c) 2008 Gabriel Peyre

options.null = 0;
rho = getoptions(options, 'rho',.8);
qsub = getoptions(options, 'qsub',2);
name = getoptions(options, 'name','');

rand('state',123456);
switch type
    case 'rand'
        sel = randperm(n^2);
        mask = zeros(n);
        mask(sel(1:round(rho*end)))=1;
    case 'mask'
        n0 = getoptions(options, 'n0',1,1);
        mask = load_image([name '-mask'], n0);
        mask = mean(mask,3);
        mask = rescale( crop(mask,n), 0,1 );
        mask = rescale(mask)<.5;
    case 'checkboard'
        x = floor(0:1/(qsub-1):n); x= x(1:n);
        [Y,X] = meshgrid(x,x);
        mask = mod( X+Y,2 )==0;
    case 'superresol'
        [Y,X] = meshgrid(0:n-1,0:n-1);
        mask = (mod(X,qsub)==0) & (mod(Y,qsub)==0);
        mask = 1-mask; 
    case {'grids' 'grids-inv'}
        [Y,X] = meshgrid(0:n-1,0:n-1);
        ngrid = getoptions(options, 'ngrid', 6);
        width = getoptions(options, 'width', 8);
        delta = n/ngrid;
        s = round(delta/2:delta:n);
        mask = zeros(n);
        for i=1:length(s);
            mask( abs(X-s(i))<=width/2 ) = 1;
        end
        mask = (mask + mask')>0;
        if strcmp(type, 'grids-inv')
            mask = 1-mask;
        end
    otherwise
        error('Unknown mask');
end