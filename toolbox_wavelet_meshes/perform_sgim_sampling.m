function gim = perform_sgim_sampling(signal_original, pos_sphere, face, n, sampling_type)

% perform_sgim_sampling - generate a geometry image from a spherical parameterization
%
%   gim = perform_sgim_sampling(signal_original, pos_sphere, face, n);
%
%   'signal_original' are typically vertex location of the original mesh, but it can 
%   also be any information you want to re-sample on a regular grid (normal, etc).
%   'pos_sphere' are location on the sphere for these points.
%   'face' is the face data structure of the mesh.
%   'n' is width of the GIM.
%
%   Uses the spherical geometry image datastructure introduced in
%       Spherical parametrization and remeshing.
%       E. Praun, H. Hoppe.
%       SIGGRAPH 2003, 340-349.
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<4
	n = 64;
end

if nargin<5
    sampling_type = 'area';
end

if size(signal_original,2)~=size(pos_sphere,2)
    error('Original and spherical meshes must be of same size.');
end

% compute sampling location on the image
disp('Computing planar sampling locations.');
posw = perform_spherial_planar_sampling(pos_sphere, sampling_type);

% perform 4-fold symmetry 
sym = { [0,1], [0,-1], [1, 0], [-1, 0] };
for i=1:length(sym)
    c = sym{i};
    if c(1)==0
        I = find(posw(2,:)*sign(c(2))>=0);
    else
        I = find(posw(1,:)*sign(c(1))>=0);
    end
    posi = posw(:,I);
    posi = perform_symmetry(posi,c);
    posw = [posw, posi];
    signal_original = [signal_original, signal_original(:,I)];
end

% crop a bit to speed up ..
m = 1.5;
I = find( posw(1,:)<m & posw(1,:)>-m & posw(2,:)<m & posw(2,:)>-m);
posw = posw(:,I);
signal_original = signal_original(:,I);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTERPOLATE THE CENTER (orignal triangulation)
% remove the faces that cross the boundary
edges = [face(1:2,:), face(2:3,:), face([1 3],:)];
pos = pos_sphere;
J = find( pos(3,edges(1,:))<=0 & pos(3,edges(2,:))<=0 &     ...
          ( pos(1,edges(1,:)).*pos(1,edges(2,:))<=0   |     ...
            pos(2,edges(1,:)).*pos(2,edges(2,:))<=0 )  );
% find the corresponding faces
J = unique( mod(J-1,size(face,2))+1 );
% find the complementary
x = zeros(size(face,2),1); x(J) = 1; I = find(x==0);
% retrive face number
face1 = face; face1 = face1(:,I);
% compute the interpolation on this sub-set using original triangulation
posn = (n-1)*(posw+1)/2+1;      % sampling location in [1,n]²
gim1 = zeros( n, n, size(signal_original,1) );
fprintf('Griding using original triangulation ');
for i=1:size(signal_original, 1)
    fprintf('.');
    gim1(:,:,i) = griddata_arbitrary( face1, posn, signal_original(i,:), n );
end


% remove doublons
if 1
[tmp,I] = unique(posw(1,:)+pi*posw(2,:));
posw = posw(:,I);
signal_original = signal_original(:,I);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTERPOLATE THE BOUNDARY (delaunay triangulation)
% sampling location in [-1,1]
x = -1:2/(n-1):1;
[Y,X] = meshgrid(x,x);
% interpolate location
gim2 = zeros( n, n, size(signal_original,1) );
fprintf('\nGriding using Delaunay triangulation ');
for i=1:size(signal_original, 1)
    fprintf('.');
    gim2(:,:,i) = griddata( posw(1,:), posw(2,:), signal_original(i,:), X, Y );
end
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MIX THE TWO
gim = gim1;
I = find( isnan(gim) );
gim(I) = gim2(I);

I = find( isnan(gim) );
gim(I) = 0;

% to keep delaunay uncomment this
% gim = gim2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% translate/scale to fit in a box [0,1]^3
for i=1:size(signal_original, 1)
    x = gim(:,:,i); x = x(:);
    m = min( x );
    gim(:,:,i) = gim(:,:,i) - m;
end
gim = rescale(gim);

function y = perform_symmetry(x,c);
% y = 2*c - x
y(1,:) = 2*c(1) - x(1,:);
y(2,:) = 2*c(2) - x(2,:);