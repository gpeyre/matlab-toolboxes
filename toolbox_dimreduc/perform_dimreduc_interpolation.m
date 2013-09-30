function A = perform_dimreduc_interpolation(pos,xy,X,options)


% full dimension
d = size(X,1);
% reduced dimension
d1= size(xy,1);
% number of points
p = size(X,2);

if size(pos,2)>1
    % animation along a path
    A = zeros( d,size(pos,2) );
    for i=1:size(pos,2)
        A(:,i) = perform_dimreduc_interpolation(pos(:,i),xy,X,options);
    end
    return;
end

% perform interpolation
options.null = 0;
if isfield(options, 'interp_mode')
    interp_mode = options.interp_mode;
else
    interp_mode = 'gaussian';
end

switch lower(interp_mode)
    case {'gaussian','recipr','nearest'}
        % estimate sigma
        sigma = sqrt( sum( std( xy' ).^2 ) );
        sigma = 4 * sigma / p.^(1/d1);
        % compute the weights
        pos = repmat( pos, [1 p] );
        dist = sqrt( sum( (xy-pos).^2,1 ) );
        if strcmp(interp_mode,'gaussian')
            w = exp( -dist.^2 / (2*sigma^2) );
        elseif strcmp(interp_mode,'nearest')
            [tmp,I] = min(dist);
            w = dist*0; w(I) = 1;
        else
            if isfield(options, 'exponent')
                eta  = options.exponent;
            else
                eta = 3;
            end
            w = 1 ./ ( (dist/sigma).^eta + 1 );
        end   
        w = w / sum(w);
        w = repmat( w, [d,1] );
        A = sum( w .* X, 2 );
        
    case 'delaunay'
        % works only in 2D
        if d1>2
            warning('Works only in 2D, ignoring other dimensions.');
        end
end
