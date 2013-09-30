function [vertex,face] = mesh_previewer(file)

% off_previewer - load a file in OFF, PLY, WRL
%   format and display the model.
%
%   [vertex,face] = mesh_previewer(file);   % retrieve faces and vertices
%   mesh_previewer(file);
%   mesh_previewer;     % open a dialog box
%
%   Use it with no argument if you want
%   to open a file browser.
%
%   See also: plot_mesh.
%
%   Copyright (c) 2004 Gabriel Peyré

h = figure;
if nargin==0
    [f, pathname] = uigetfile({'*.off;*.ply;*.wrl;*.smf;*.png;*.jpg;*.gim','*.off,*.ply,*.wrl,*.smf,*.png,*.png,*.gim Files'},'Pick a file');
    file = [pathname,f];
end

ext = file(end-2:end);
ext = lower(ext);

if strcmp(ext,'gim') || strcmp(ext,'png') || strcmp(ext,'jpg')
    gim = read_gim(file);
    plot_geometry_image(gim);
    if nargin>0
        [vertex,face] = convert_gim2mesh(gim);
    end
else
    [vertex,face] = read_mesh(file);
    plot_mesh(vertex,face);
end
shading interp;
cameramenu;