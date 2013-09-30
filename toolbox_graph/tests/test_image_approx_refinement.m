% test for adaptive triangulation using refinement

n = 256;
name = 'disk';
name = 'triangle';
name = 'lena';
M = load_image(name,n);
M = rescale(M);

[Y,X] = meshgrid([1 n], [1 n]);
vertex = cat(1, X(:)', Y(:)');
vertex = vertex+randn(size(vertex))*.001;
% vertex(:,end+1) = [1; 10];
% vertex(:,end+1:end+5) = rand(2,5)*(n-1)+1;

[Xi,Yi] = meshgrid(1:n,1:n);

niter = 1000;

for it=1:niter
    
    eta = 1e-3;

    % compute voronoi regions
    p = size(vertex,2);
    Q = griddata(vertex(1,:)+randn(1,p)*eta, vertex(2,:)+randn(1,p)*eta, 1:p, Xi,Yi, 'nearest');
    % face = compute_voronoi_triangulation(Q, vertex);
    

    % extract steiner points
    Q1 = zeros(n+2)-1;
    Q1(2:end-1,2:end-1) = Q;
    V = [];
    v = Q1(1:end-1,1:end-1); V = [V v(:)];
    v = Q1(2:end,1:end-1); V = [V v(:)];
    v = Q1(1:end-1,2:end); V = [V v(:)];
    v = Q1(2:end,2:end); V = [V v(:)];
    V = sort(V,2);
    d = (V(:,1)~=V(:,2)) + (V(:,2)~=V(:,3)) + (V(:,3)~=V(:,4));
    I = find(d>=2);
    [vx,vy] = ind2sub([n+1 n+1], I);
    vertexv = cat(1,vy',vx');
    
    if 0
    hold on;
    imageplot(L); 
    plot(vertex(1,:), vertex(2,:), '.');
    plot(vertexv(1,:), vertexv(2,:), 'r.');
    axis tight;
    hold off;
    end

    % compute voronoi cells associated to voronoi points
    q = size(vertexv,2);
    warning off;
    L = griddata(vertexv(1,:)+randn(1,q)*eta, vertexv(2,:)+randn(1,q)*eta, 1:q, Xi,Yi, 'nearest');
    warning on;

    % compute approximation error
    u = unique(L(:));
    nu = size(vertexv,2);
    e = zeros(nu,1);
    for k=1:nu
        v = M(L==k);
        if not(isempty(v))
            mu = mean(v);
            e(k) = sum( (v-mu).^2 );
        end
    end
    % add point
    [tmp,i] = max(e);
    vertex = cat(2, vertex, vertexv(:,i));
    
    face = delaunay(vertex(1,:),vertex(2,:))';    
    clf;
    hold on;
    imageplot(M);% colormap jet;
    options.col = 'b';
    plot_graph( triangulation2adjacency(face), vertex, options );
    h = plot(vertex(1,:), vertex(2,:), 'r.');
    set(h, 'MarkerSize', 20);
    hold off;
    axis ij;
    drawnow;
end
