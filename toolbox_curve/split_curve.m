function c_list1 = split_curve(c_list, split_strength)

% split_curve - split where high curvature varitation occurs.
%
%   c_list1 = split_curve(c_list, split_strength);
%
%   Copyright (c) 2004 Gabriel Peyré

if ~iscell(c_list)
    c_list = {c_list};
end
if nargin<2
    split_strength = 1/1.5;
end

% segment curves
nb = length(c_list);
c_list1 = {};
j = 0;
for i=1:nb
    c = c_list{i};
    % smooth the curve
    x = c(1,:);
    y = c(2,:);
    h = compute_gaussian_filter(5,0.2); % h = [1/4,1/2,1/4];
    x = perform_convolution(x,h);
    y = perform_convolution(y,h);
    ch = [x';y']; 
    % compute curvature
    dx = compute_diff(x);
    dy = compute_diff(y);
    s = sqrt( dx.^2 + dy.^2 );
    dx = dx./s; dy = dy./s; 
    % compute dot product
    p = length(dx);
    curv = zeros(p,1);
    for k=1:p
        if k>1 && k<p
            curv(k) = dx(k-1)*dx(k+1)+dy(k-1)*dy(k+1);
        elseif k==1
            curv(k) = dx(k)*dx(k+1)+dy(k)*dy(k+1);
        elseif k==p
            curv(k) = dx(k-1)*dx(k)+dy(k-1)*dy(k);            
        end
    end
    % extract component with low curvature variation
    curv = 1./abs(curv);
    f = curv<mean(curv)/split_strength;
    f = cumsum( abs(diff(f)) ); % one step by constant compoenent
    for t=0:max(f)
        I = find(f==t);
        if length(I)>3
            j = j+1;
            c_list1{j} = c(:,I);
        end
    end
end