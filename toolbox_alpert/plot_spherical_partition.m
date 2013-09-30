function plot_spherical_partition(xyz, part)

% plot_spherical - plot clustered points on a sphere.
%
%   plot_spherical_partition(xyz, part);
%
%   Copyright (c) 2004 Gabriel Peyré

a = xyz(1,:);
b = xyz(2,:);    
c = xyz(3,:);

% plot a nice sphere
p = 20;
t = 0:1/(p-1):1;
[th,ph] = meshgrid( t*pi,t*2*pi );
x = cos(th);
y = sin(th).*cos(ph);
z = sin(th).*sin(ph);

hold on;
surf(x,y,z, z.*0+1);

shading interp;
lighting gouraud;
camlight infinite; 
camproj('perspective');
axis square; 
axis off;

for i=1:length(part)
    I = part{i};
    str = [get_color_from_index(i), '.'];
    plot3(a(I),b(I),c(I), str);
end

axis square;
hold off;