function plot_spherical_triangulation(svertex,face, options)


% plot_spherical_triangulation - display a nice spherical triangulation
%
%   plot_spherical_triangulation(svertex,face,options);
%
%   Copyright (c) 2008 Gabriel Peyre

options.null = 0;


hold on;
% draw a background sphere
[X,Y,Z] = sphere(30); 
surf(X,Y,Z);
shading interp; lighting gouraud
camlight;

f = getoptions(options, 'target_face', []);
f = face(:,f);

edges = compute_edges(face);
ne = size(edges,2);
npoints = getoptions(options, 'npoints', 30); % number of point on each arc
lw = getoptions(options, 'lw', 3);
for k=1:ne
    i = edges(1,k); j = edges(2,k);
    a = svertex(:,i); b = svertex(:,j);
    if not(ismember(i,f)) || not(ismember(j,f))
        plot_arc(a, b, npoints, lw, 'b');
    end
end
axis tight; axis equal; axis off;

rho = 1.02;
if not(isempty(f))
    for i=1:3
        j = mod(i,3)+1;
        c = svertex(:,f(i)); d = svertex(:,f(j));
        h = plot3(c(1)*rho, c(2)*rho, c(3)*rho, '.g'); 
        set(h, 'MarkerSize', 40);
        plot_arc(c, d, npoints, lw+1, 'r');
    end
end

%%
function plot_arc(a, b, npoints, lw, col)

t = repmat(linspace(0,1,npoints), [3 1]);

% arc curve
c = t.*repmat(a, [1 npoints]) + (1-t).*repmat(b, [1 npoints]);
% project on sphere
d = sqrt(sum(c.^2)); d(d<1e-10) = 1;
c = 1.005 * c./repmat(d, [3 1]);
h = plot3(c(1,:), c(2,:), c(3,:), col);
set(h, 'LineWidth', lw);