function theta = compute_vf_angle(v,dir)

if nargin<2
    dir = 1;
    if size(v,3)==1
        dir=-1;
    end
end
if dir==1
    theta = atan2(v(:,:,2),v(:,:,1));
else
    theta = cat(3,cos(v),sin(v));
end