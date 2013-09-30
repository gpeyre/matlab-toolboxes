function plot_flattened_dataset(xy,M,k);

% plot_flattened_dataset - nice 2D layout of image according to sampling location xy
%
% plot_flattened_dataset(xy,M,k);
%
%	M is a matrix of size (a x b x nbr_images), a collection of images of size (a,b)
%	xy should be of size (nbr_sample x 2)
%	k is the number of images stiched in each direction.
%
%   Copyright (c) 2005 Gabriel Peyré

if nargin<3
	k = 30;
end

if size(xy,2)>size(xy,1)
    xy = xy';
end


a = size(M,1);
b = size(M,2);
n = size(M,3);


e = 1/(2*k);
% plot result
for i=1:2
    xy(:,i) = rescale(xy(:,i), e, 1-e );
end

A = zeros(k*b,k*a) + 1; % mmax(M(:));
for x=1:k
    for y=1:k
        selx = (x-1)*b+1:x*b;
        sely = (y-1)*a+1:y*a;
        % center point
        cx = e + (x-1)/k;
        cy = e + (y-1)/k;
        % distance to center
        d = max( abs( xy(:,1)-cx ), abs( xy(:,2)-cy ) );
        % find closest point
        [v,I] = min(d);
        if v<=e
            A(selx,sely) = rescale( M(end:-1:1,:,I)' );
        end
    end
end


n = size(A,1);
q = size(A,2);
xy(:,1) = xy(:,1)*n;
xy(:,2) = xy(:,2)*q;

plot_lines = 1;

hold on;
imagesc( [0,n-1]+0.5,[0,q-1]+0.5, A' );
colormap gray(256);
axis tight; axis square; axis off;
scatter(xy(:,1),xy(:,2),12,'r','filled');
axis xy;
if plot_lines
    for x=0:b:n
        line([x x], [0 q]);
    end
    for x=0:a:q
        line([0 n], [x x]);
    end
end
hold off;