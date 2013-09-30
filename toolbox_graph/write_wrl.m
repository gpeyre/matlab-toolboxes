function write_wrl(filename, vertex, face);

% write_wrl - write a mesh to a VRML file
%
%   write_wrl(filename, vertex, face);
%
%   vertex must be of size [n,3]
%   face must be of size [p,3]
%
%   Copyright (c) 2003 Gabriel Peyré


if size(vertex,2)~=3
    vertex = vertex';
end
if size(vertex,2)~=3
    error('vertex is not of correct format.');
end

if size(face,2)~=3
    face = face';
end
if size(face,2)~=3
    error('face is not of correct format.');
end

fid = fopen(filename,'wt');
if( fid==-1 )
    error('Can''t open the file.');
    return;
end

% header
fprintf(fid, '#VRML V2.0 utf8\nShape { \nappearance Appearance { \nmaterial Material { \n}\n} \ngeometry IndexedFaceSet { \ncoord Coordinate { \npoint [\n');

% write the points
for i=1:size(vertex,1)
    pos = vertex(i,:);
    fprintf(fid, '%f %f %f,\n', pos(1), pos(2), pos(3));
end

fprintf(fid, '] \n} \ncoordIndex [ \n');

% write the faces
for i=1:size(face,1)
    f = face(i,:)-1;    % index start at 0
    fprintf(fid, '%d %d %d -1,\n', f(1), f(2), f(3));
end

fprintf(fid, '] \n} \n}\n');

fclose(fid);