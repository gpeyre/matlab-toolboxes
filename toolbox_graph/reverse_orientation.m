function face1 = reverse_orientation(face)

% reverse_orientation - reverse the orientation
%   of the face of a mesh.
%
%   face1 = reverse_orientation(face);
%
%   Copyright (c) 2003 Gabriel Peyré

nface = size(face);
face1 = face;

for i=1:nface
    f = face1(i,:);
    face1(i,:) = f(3:-1:1);
end
