function face = perform_triangle_flipping(face, flips, options)

% perform_triangle_flipping - apply a sequence of flips to a triangulation
%
%   face = perform_triangle_flipping(face, flips, options);
%
%   Flip that are not topologically valid (either on a boundary edge or a
%   flip that creates a double face) are skipped.
%
%   flips(:,k) is the edge [i;j] that is flipped at step k.
%
%   You can provide your own listing of edge via options.edge.
%
%   Copyright (c) 2008 Gabriel Peyre

options.null = 0;
warning_boundary = getoptions(options, 'warning_boundary', 1);

% edge = compute_edges(face);

% face attached to edges
e2f = compute_edge_face_ring(face);

for k=1:size(flips,2)
    % edge to flip
    e = flips(:,k); % edge(:,flips(k));
    i = e(1); j = e(2);
    if i==j
        warning('Problem');
    end
    % adjacent edges
    f1 = e2f(i,j); f2 = e2f(j,i);
    if f1<=0 || f2<=0
        if warning_boundary
            warning('Cannot flip boundary edges');
        end
    else
        % find the two other points
        k1 = find_other( face(:,f1), i,j );
        k2 = find_other( face(:,f2), i,j );
        % new faces
        t1 = [k1;k2;i];
        t2 = [k1;k2;j];
        % check wether the faces are not yet in the triangulation
        I = intersect(sort(face)',sort([t1 t2])','rows');
        if isempty(I)
            face(:,f1) = t1;
            face(:,f2) = t2;
            % update e2f
            % todo : speed up this ...
            e2f = compute_edge_face_ring(face);
        else
            if warning_boundary
                warning('Invalid flip (create already existing face)');
            end
        end
    end
end

%%%%%%%%%
function f = find_other( f, i,j )

I = find(f==i);
if length(I)~=1
    error('Problem');
end
f(I) = [];
I = find(f==j);
if length(I)~=1
    error('Problem');
end
f(I) = [];