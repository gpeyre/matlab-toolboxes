function X = perform_lloyd_iteration(X,options)

% perform_lloyd_iteration - perform lloyd smoothing
%
%   X = perform_lloyd_iteration(X,options);
%
%   X is a (2,n) matrix of n points in [0,1]x[0,1]
%   The algorithm perform a step of Lloy relaxation to evenly distribute
%   the points in the square.
%
%   You can give a local density using an image options.density
%
%   Copyright (c) Gabriel Peyr   

options.null = 0;
niter = getoptions(options, 'niter',1);
D = getoptions(options, 'density',[]);

n = size(X,2);
ndensity = size(D,1);

%% add point at infinity
a = 5;
X(:,end+1) = [a;a];
X(:,end+1) = [-a;a];
X(:,end+1) = [-a;-a];
X(:,end+1) = [a;-a];

for k=1:niter

    [V,C] = voronoin(X');
    if isempty(D)
        for j=1:n
            X(:,j) = compute_polygon_center( V(C{j},:)' );
        end
    else
        % create a triangulation
        T = []; v = [];
        nc = size(V,1);
        vertex = [V' X];
        for j=1:n
            vertex(:,j+nc) = mean( V(C{j},:), 1)';
            qc = length(C{j});
            c = [C{j} C{j}(1)];
            T(:,end+1:end+qc) = [c(1:end-1); c(2:end); repmat(j+nc, 1, qc)];
            v(end+1:end+qc) = j;
        end
        % index map
        vertex = vertex*(ndensity-1)+1;
        options.verb = 0;
        M = griddata_arbitrary(T,vertex,v,ndensity, options);
        % compute center of gravity
        x = linspace(0,1,ndensity);
        [Yp,Xp] = meshgrid(x,x);
        for j=1:n
            I = find(M==j);
            if not(isempty(I))
                d = D(I); % density weights
                d = rescale(d,.1,1);
                X(:,j) = [ sum(Xp(I).*d); sum(Yp(I).*d) ] / sum(d);
            else
                X(:,j) = [1/2 1/2];
            end
        end
    end
    X(:,1:n) = clamp(X(:,1:n),0,1);

end

X = X(:,1:n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = compute_polygon_area(P)

P(:,end+1) = P(:,1);
x = P(1,:);
y = P(2,:);

A = 1/2 * sum( x(1:end-1).*y(2:end)-y(1:end-1).*x(2:end) );
A = abs(A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = compute_polygon_center(P)

P = clamp(P,0,1);
A = compute_polygon_area(P);
P(:,end+1) = P(:,1);
x = P(1,:);
y = P(2,:);
c(1) = sum( (x(1:end-1)+x(2:end)) .* ( x(1:end-1).*y(2:end)-y(1:end-1).*x(2:end) ) );
c(2) = sum( (y(1:end-1)+y(2:end)) .* ( x(1:end-1).*y(2:end)-y(1:end-1).*x(2:end) ) );
if A>0
    c = c/(6*A);
else
    c = [1/2 1/2];
end
c = abs(c);