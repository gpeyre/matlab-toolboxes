function plot_spherical_function(a,b,c,v)

% plot_spherical_function - plot function defined on a sphere.
%
%   plot_spherical_function(a,b,c,[v]);
% or
%   plot_spherical_function(points,[v]);
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin==1 || nargin==2
    if nargin==1
        v = [];
    else
        v = b;
    end
    b = a(2,:);    
    c = a(3,:);
    a = a(1,:);
end
if nargin==3
    v = [];
end

v = v(:);

% prevision for the drawing
p = 50;

hold on;
if 0
    [x,y,z] = sphere(p);
else
    t = 0:1/(p-1):1;
    [th,ph] = meshgrid( (t-0.5)*pi,(t-0.5)*2*pi );
    x = sin(th);
    y = cos(th).*sin(ph);
    z = cos(th).*cos(ph);
end

if ~isempty(v)
    
    [Ph,Th,R] = cart2sph(a,b,c);    
    Ph = [Ph,Ph+2*pi,Ph-2*pi];
    Th = [Th,Th,Th];
    v1 = [v',v',v'];
    
    % use scattered data interpolation
    V = griddata(Th,Ph,v1,th,ph);
    V = rescale(V,-1,1);
    
    epsi = 0.2;
    if 0
        x = x + epsi*V;
        y = z + epsi*V;
        z = z + epsi*V;
        
        a = a + epsi*v;
        b = b + epsi*v;
        c = c + epsi*v;
    end
end

colormap jet(256);
if exist('V')
    surf(x,y,z, rescale(V));
else
    surf(x,y,z, z.*0+1);
end

shading interp;
lighting gouraud;
camlight infinite; 
camproj('perspective');
axis square; 
axis off;

plot3(a,b,c, '.');

axis square;
hold off;