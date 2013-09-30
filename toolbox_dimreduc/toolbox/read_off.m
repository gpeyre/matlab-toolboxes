function [vertex,face] = read_off(filename, verbose)

% read_off - read data from OFF file.
%
%   [vertex,face] = read_off(filename, verbose);
%
%   'vertex' is a 'nb.vert x 3' array specifying the position of the vertices.
%   'face' is a 'nb.face x 3' array specifying the connectivity of the mesh.
%
%   Copyright (c) 2003 Gabriel Peyré

if nargin<2
    verbose = 1;
end

fid = fopen(filename,'r');
if( fid==-1 )
    error('Can''t open the file.');
    return;
end

str = fgets(fid);   % -1 if eof
if ~strcmp(str(1:3), 'OFF')
    error('The file is not a valid OFF one.');    
end

str = fgets(fid);
[a,str] = strtok(str); nvert = str2num(a);
[a,str] = strtok(str); nface = str2num(a);



[A,cnt] = fscanf(fid,'%f %f %f', 3*nvert);
if cnt~=3*nvert
    warning('Problem in reading vertices.');
end
A = reshape(A, 3, cnt/3);
vertex = A;
% read Face 1  1088 480 1022
[A,cnt] = fscanf(fid,'%d %d %d %d\n', 4*nface);
if cnt~=4*nface
    warning('Problem in reading faces.');
end
A = reshape(A, 4, cnt/4);
face = A(2:4,:)+1;


fclose(fid);





return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% old code

str = ['Load mesh ', filename, ' : ', num2str(nvert), ' vertices, ', num2str(nface), ' faces.'];
str = strrep(str,'_','\_'); % problem when displaying '_' (LaTeX special char)
str = strrep(str,'\','/'); % problem when displaying '\' (LaTeX special char)
if verbose
    h = waitbar(0,str);
end

vertex = zeros(nvert,3);
face = zeros(nface,3);

% loading vertex
i = 0;
while i<nvert
    if verbose
        waitbar(i/(nvert+nface));
    end
    str = fgets(fid);
    if ~strcmp(strtrim(str), '')
        i = i+1;
        [a,str] = strtok(str); x = str2num(a);
        [a,str] = strtok(str); y = str2num(a);
        [a,str] = strtok(str); z = str2num(a);
        vertex(i,:) = [x y z];
    end
end

% loading face
for i=1:nface
    if verbose
        waitbar((i+nvert)/(nvert+nface));
    end
    str = fgets(fid);
    [a,str] = strtok(str); nb = str2num(a);
    if  nb~=3
        error('The file is not a valid OFF one.');
    end
    [a,str] = strtok(str); x = str2num(a);
    [a,str] = strtok(str); y = str2num(a);
    [a,str] = strtok(str); z = str2num(a);
    face(i,:) = [x y z]+1;
end
if verbose
    close(h);
end

fclose(fid);