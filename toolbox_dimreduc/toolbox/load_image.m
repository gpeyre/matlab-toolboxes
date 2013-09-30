function M = load_image(type, n, options)

% load_image - load benchmark images.
%
%   M = load_image(name, n, options);
%
%   name can be:
%   Synthetic images:
%       'chessboard1', 'chessboard', 'square', 'disk', 'quaterdisk', '3contours', 'line',
%       'line_vertical', 'line_horizontal', 'line_diagonal', 'line_circle',
%       'parabola', 'sin', 'phantom'
%   Natural images:
%       'boat', 'lena', 'goldhill', 'mandrill', 'maurice', 'polygons_blurred'.
%   
%   Copyright (c) 2004 Gabriel Peyré

if nargin<2
    n = 512;
end
options.null = 0;

type = lower(type);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters for geometric objects
eta = 0.1;              % translation
gamma = 1/sqrt(2);      % slope
if isfield( options, 'eta' )
    eta = options.eta;
end
if isfield( options, 'gamma' )
    eta = options.gamma;
end
if isfield( options, 'radius' )
    radius = options.radius;
end
if isfield( options, 'center' )
    center = options.center;
end
if isfield( options, 'center1' )
    center1 = options.center1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for the line, can be vertical / horizontal / diagonal / any
if strcmp(type, 'line_vertical')
    eta = 0.5;              % translation
    gamma = 0;      % slope
elseif strcmp(type, 'line_horizontal')
    eta = 0.5;              % translation
    gamma = Inf;      % slope
elseif strcmp(type, 'line_diagonal')
    eta = 0;              % translation
    gamma = 1;      % slope
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for some blurring
sigma = 0;
if isfield(options, 'sigma')
    sigma = options.sigma;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch type
    case 'chessboard1'
        x = -1:2/(n-1):1;
        [Y,X] = meshgrid(x,x);
        M = (2*(Y>=0)-1).*(2*(X>=0)-1);
        
    case 'chessboard'
        if ~isfield( options, 'width' )
            width = round(n/16);
        else
            width = options.width;
        end
        [Y,X] = meshgrid(0:n-1,0:n-1);
        M = mod( floor(X/width)+floor(Y/width), 2 ) == 0;
        
    case 'square'
        if ~isfield( options, 'radius' )
            radius = 0.6;
        end
        x = -1:2/(n-1):1;
        [Y,X] = meshgrid(x,x);
        M = max( abs(X),abs(Y) )<radius;
        
    case {'line', 'line_vertical', 'line_horizontal', 'line_diagonal'}
        x = 0:1/(n-1):1;
        [Y,X] = meshgrid(x,x);
        if gamma~=Inf
            M = (X-eta) - gamma*Y < 0;
        else
            M = (Y-eta) < 0;
        end
        
    case 'disk'
        if ~isfield( options, 'radius' )
            radius = 0.35;
        end
        if ~isfield( options, 'center' )
            center = [0.5, 0.5];    % center of the circle
        end
        x = 0:1/(n-1):1;
        [Y,X] = meshgrid(x,x);
        M = (X-center(1)).^2 + (Y-center(2)).^2 < radius^2;
        
    case 'quarterdisk'
        if ~isfield( options, 'radius' )
            radius = 0.95;
        end
        if ~isfield( options, 'center' )
            center = -[0.1, 0.1];    % center of the circle
        end
        x = 0:1/(n-1):1;
        [Y,X] = meshgrid(x,x);
        M = (X-center(1)).^2 + (Y-center(2)).^2 < radius^2;
        
    case 'fading_contour'
        if ~isfield( options, 'radius' )
            radius = 0.95;
        end
        if ~isfield( options, 'center' )
            center = -[0.1, 0.1];    % center of the circle
        end
        x = 0:1/(n-1):1;
        [Y,X] = meshgrid(x,x);
        M = (X-center(1)).^2 + (Y-center(2)).^2 < radius^2;
        theta = 2/pi*atan2(Y,X);
        h = 0.5;
        M = exp(-(1-theta).^2/h^2).*M;
        
    case '3contours'        
        radius = 1.3;
        center = [-1, 1];
        radius1 = 0.8;
        center1 = [0, 0];
        x = 0:1/(n-1):1;
        [Y,X] = meshgrid(x,x);
        f1 = (X-center(1)).^2 + (Y-center(2)).^2 < radius^2;
        f2 = (X-center1(1)).^2 + (Y-center1(2)).^2 < radius1^2;
        M = f1 + 0.5*f2.*(1-f1);
        
    case 'line_circle'

        gamma = 1/sqrt(2);

        x = linspace(-1,1,n);
        [Y,X] = meshgrid(x,x);
        M1 = double( X>gamma*Y+0.25 );
        M2 = X.^2 + Y.^2 < 0.6^2;
        M = 20 + max(0.5*M1,M2) * 216;
    
    case {'bwhorizontal','bwvertical','bwcircle'}
        
        [Y,X] = meshgrid(0:n-1,0:n-1);
        if strcmp(type, 'bwhorizontal')
            d = X;
        elseif strcmp(type, 'bwvertical')
            d = Y;
        elseif strcmp(type, 'bwcircle')
            d = sqrt( (X-(n-1)/2).^2 + (Y-(n-1)/2).^2 );
        end
        if isfield(options, 'stripe_width')
            stripe_width = options.stripe_width;
        else
            stripe_width = 5;
        end
        if isfield(options, 'black_prop')
            black_prop = options.black_prop;
        else
            black_prop = 0.5;
        end
        M = double( mod( d/(2*stripe_width),1 )>=black_prop );
        
    case 'parabola'
        
        % curvature
        if isfield(options, 'c')
            c = options.c;
        else
            c = 0.1;
        end
        % angle
        if isfield(options, 'theta');
            theta = options.theta;
        else
            theta = pi/sqrt(2);
        end
        x = -0.5:1/(n-1):0.5;
        [Y,X] = meshgrid(x,x);
        Xs = X*cos(theta) + Y*sin(theta);
        Y =-X*sin(theta) + Y*cos(theta); X = Xs;
        M = Y>c*X.^2; 
        
    case 'sin'
        
        [Y,X] = meshgrid(-1:2/(n-1):1, -1:2/(n-1):1);
        M = Y >= 0.6*cos(pi*X);
        M = double(M);

    case 'phantom'
        
        M = phantom(n);
        
    case 'periodic_bumps'
        
        if isfield(options, 'nbr_periods')
            nbr_periods = options.nbr_periods;
        else
            nbr_periods = 8;
        end
        if isfield(options, 'theta')
            theta = options.theta;
        else
            theta = 1/sqrt(2);
        end
        if isfield(options, 'skew')
            skew = options.skew;
        else
            skew = 1/sqrt(2);
        end
        
        A = [cos(theta), -sin(theta); sin(theta), cos(theta)];
        B = [1 skew; 0 1];
        T = B*A;
        x = (0:n-1)*2*pi*nbr_periods/(n-1);
        [Y,X] = meshgrid(x,x);
        pos = [X(:)'; Y(:)'];
        pos = T*pos;
        X = reshape(pos(1,:), n,n);
        Y = reshape(pos(2,:), n,n);
        M = cos(X).*sin(Y);      
        
    case 'noise'
        if isfield(options, 'sigma')
            sigma = options.sigma;
        else
            sigma = 1;
        end
        M = randn(n);
        
    otherwise
        ext = {'gif', 'png', 'jpg', 'bmp', 'tiff', 'pgm'};
        for i=1:length(ext)
            name = [type '.' ext{i}];
            if( exist(name) )
                M = imread( name );
                M = double(M);
                if n~=size(M, 1) && nargin>=2
                    M = image_resize(M,n,n);
                end
                return;
            end
        end
        error( ['Image ' type ' does not exists.'] );
end

M = double(M);

if sigma>0
    h = ones(sigma);
    h = conv2(h,h);
    h = h/sum(h(:));
    M = conv2(M, h, 'same');
end

M = rescale(M) * 256;