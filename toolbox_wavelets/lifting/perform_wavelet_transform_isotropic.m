function y = perform_wavelet_transform_isotropic(x,J,options)

% perform_wavelet_transform_isotropic - multidimensional isotropic wavelet transform (aka classical ND transform).
% 
%   y = perform_wavelet_transform_isotropic(x,J,options);
%
%   Perform a isotropic transform of the array x, on J scale.
%   'x' can be of any dimension. Recalls that the results 
%   are computed 'in place', not using Mallat's ordering.
%   Please see 'fwt_lifting' for an overview of what 
%   'in place' means in 1D.
%
%   'options' is an (optional) structure that can contain:
%       - 'verb' : control verbosity.
%       - 'type' : the type of 1D wavelet transform used (a string), e.g.
%           'perform_79_transform'.
%
%   Copyright (c) 2003 Gabriel Peyré

if J<1
    y = x;
    return; 
end

dim = size(x);
d = length(dim);
n = prod(dim);

% remove empty dimension
while( dim(d)==1 )
    d = d-1;
    dim = dim(1:d);    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3
    options.null = 0;   % force creation of options
end

if isfield(options,'type')   % kind of transform
    type = options.type;
else
    type = 'perform_79_transform';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begining of the code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2
    J = floor(log2(min(dim)));
end

y = x;

% perform 1 step of transform for each dimension
for s=1:d
    
    if d~=1
        curd = dim([1:(s-1),s+1:d]);
    else
        curd = [1];
    end
    
    for t=1:prod(curd)
    
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % index computations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if d~=1
            % turn t into an index
            [out{1:d-1}] = ind2sub(curd,t);
            ind_iso = cell2mat(out);
            % put the vector in v
            v = array_get_dimension(y, ind_iso, s);
        else
            ind_iso = 1;
            v = y;
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % perform transform
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % perform the transform on 1 scale only
        str = sprintf('v = %s(v,1);', type);
        eval(str);        
        % assign the result
        if d~=1
            y = array_set_dimension( y,ind_iso,s,v );
        else
            y = v;
        end
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% launch the recursion on a sub square
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if J>1
    
    % selection string for coarse channel
    sel = [];
    for i=1:d
        sel = [sel, sprintf('1:2:%d', dim(i))];
        if i~=d
            sel = [sel, ','];    
        end
    end
    
    % perform extraction of coarse channels
    str = ['M = y(', sel, ');'];
    eval(str);
    
    % set up options for recursion, mainly discontinuity extraction
    clear options; options.null = 0; % force creation
    options.type = type;
    
    % perform transform
    M = perform_wavelet_transform_isotropic(M,J-1,options);
    
    % perform assignement
    str = ['y(', sel, ') = M;'];
    eval(str);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
