function [v,f] = mesh_flip_faces(arg1, arg2, arg3)

if nargin==1
    name = arg1;
    [v,f] = read_mesh(name);
else
    v = arg1;
    f = arg2
end

if size(f,1)~=3
    f = f';
end

f = [f(1,:); f(3,:); f(2,:)];

if nargin==1
    write_mesh(name, v,f);
end