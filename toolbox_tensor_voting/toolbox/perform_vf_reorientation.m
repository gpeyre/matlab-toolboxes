function w = perform_vf_reorientation(v, method, w)

% perform_vf_reorientation - try to reorient the vf.
%
%   w = perform_vf_reorientation(v, method);
%
%   'method' can be 'xproj', 'yproj', 'circproj' or 'localproj' or
%   'custproj'.
%
%   For 'custproj' you have to provide an additional vector w.
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<2
    method = 'xproj';
end


if strcmp(method, 'localproj')
    % special case.
    n = size(v,1);
    p = size(v,2);
    for i=1:n
        for j=1:n
            m = zeros(1,1,2);
            if i>1
                m = m + v(i-1,j,:);
            end
            if j>1
                m = m + v(i,j-1,:);
            end
            if i>1 && j>1
                m = m + v(i-1,j-1,:);
            end
            s = dot( m, v(i,j,:) );
            if s>0
                v(i,j,:) = v(i,j,:) * sign(s);
            end
                
        end
    end        
    
    w = v;
    return;
end

switch lower(method)
    case 'xproj'    
        s = v(:,:,1);
    case 'yproj'
        s = v(:,:,2);
    case 'circproj'
        n = size(v,1);
        p = size(v,2);
        [Y,X] = meshgrid(0:p-1, 0:n-1);
        s = v(:,:,1).*X + v(:,:,2).*Y;
    case 'custproj'
        s = v(:,:,1)*w(1) + v(:,:,2)*w(2);
    otherwise
        error('Unknown method');
end


w = v;
w(:,:,1) = v(:,:,1).*sign(s);
w(:,:,2) = v(:,:,2).*sign(s);
    
    