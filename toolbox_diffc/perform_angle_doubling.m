function v1 = perform_angle_doubling(v,dir)

% perform_angle_doubling - double the angle of a vector field
%
%   v1 = perform_angle_doubling(v,dir);
%
%   Copyright (c) 2008 Gabriel Peyre

theta = atan2(v(:,:,2),v(:,:,1));
if dir==1
    v1 = cat(3, cos(2*theta), sin(2*theta));
else
    v1 = cat(3, cos(theta/2), sin(theta/2));    
end