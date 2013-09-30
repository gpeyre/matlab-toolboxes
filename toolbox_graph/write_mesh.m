function write_mesh(file, vertex, face)

% write_mesh - read data to OFF, PLY, SMF or WRL file.
%
%   write_mesh(file, vertex, face);
%
%   'vertex' is a 'nb.vert x 3' array specifying the position of the vertices.
%   'face' is a 'nb.face x 3' array specifying the connectivity of the mesh.
%
%   Copyright (c) 2005 Gabriel Peyré

ext = file(end-2:end);
ext = lower(ext);
if strcmp(ext, 'off')
    write_off(file, vertex, face);
elseif strcmp(ext, 'ply')
    write_ply(file, vertex, face);
elseif strcmp(ext, 'smf')
    write_smf(file, vertex, face);
elseif strcmp(ext, 'wrl')
    write_wrl(file, vertex, face);
elseif strcmp(ext, 'obj')
    write_obj(file, vertex, face);
else
    error('Unknown extension.');    
end