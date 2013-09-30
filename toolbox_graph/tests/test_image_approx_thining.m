% test for adaptive triangulation using thining

n = 200;
name = 'disk';
name = 'lena';
name = 'triangle';
M = load_image(name,n);
M = rescale(M);


options.n = 50;
[vertex,face] = compute_base_mesh('square', 0, options);
vertex = vertex*(n-1)+1;


niter = 1000;

for i=1:niter

    m = size(face,2);
    p = size(vertex,2);

    % interpolation on each triangle
    Q = griddata_arbitrary(face,vertex,1:m,n);

    % compute error
    e = zeros(m,1);
    for i=1:m
        v = M(Q==i);
        if not(isempty(v))
            e(i) = sum( (v-mean(v)).^2 );
        end
    end
    [tmp,i] = min(e);
    
    
    % face to remove
    f = face(:,i);
    f = sort(f); f = f(end:-1:1);
    face(:,i) = [];

    % old vertices
    v = vertex(:,f);
    % new vertex
    v1 = mean(v,2);
    vertex(:,end+1) = v1;
    % remove vertices
    vertex(:,f) = [];

    if 1
        % use delaunay
        face = delaunay(vertex(1,:),vertex(2,:))';
    else
        % retriangulate by hand
        for k=1:3
            face(face==f(k)) = p+1;
        end
        for k=1:3
            face(face>f(k)) = face(face>f(k))-1;
        end
        % remove flat faces
        s = diff(sort(face));
        I = find( s(1,:)==0 | s(2,:)==0 );
        face(:,I) = [];
    end

    % display
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