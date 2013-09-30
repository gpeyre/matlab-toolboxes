function Y = perform_synthesis_quilting(X, tilesize, n, overlap, err, initial_patch)

% perform_synthesis_quilting - perform image synthesis
%
%   Y = perform_synthesis_quilting(X, tilesize, n, overlap, err, initial_patch);
%
%   The method is described in
%        ``Image Quilting for Texture Synthesis and Transfer''
%       Alexei A. Efros and William T. Freeman
%       Proceedings of SIGGRAPH '01, Los Angeles, California, August, 2001.
%
%   The implementation is by Christopher DeCoro
%   http://www.cs.princeton.edu/~cdecoro/imagequilting/
%
%   X:  The source image to be used in synthesis
%   tilesize:   the dimensions of each square tile.  Should divide size(X) evenly
%   n:  The number of tiles to be placed in the output image, in each dimension
%   overlap: The amount of overlap to allow between pixels (def: 1/6 tilesize)
%   err: used when computing list of compatible tiles (def: 0.1)

X = double(X);

if( length(size(X)) == 2 )
    X = repmat(X, [1 1 3]);
elseif( length(size(X)) ~= 3 )
    error('Input image must be 2 or 3 dimensional');
end

simple = 0;

if( nargin < 5 )
    err = 0.002;
end

if( nargin < 4 )
    overlap = round(tilesize / 6);
end
if nargin<6
    initial_patch = [];
end

use_draw = 1;
use_verb = 0;

if( overlap >= tilesize )
    error('Overlap must be less than tilesize');
end

destsize = n * tilesize - (n-1) * overlap;

Y = zeros(destsize, destsize, 3);

for i=1:n,
    for j=1:n,
        progressbar(j+(i-1)*n, n^2);
        startI = (i-1)*tilesize - (i-1) * overlap + 1;
        startJ = (j-1)*tilesize - (j-1) * overlap + 1;
        endI = startI + tilesize -1 ;
        endJ = startJ + tilesize -1;

        %Determine the distances from each tile to the overlap region
        %This will eventually be replaced with convolutions
        distances = zeros( size(X,1)-tilesize, size(X,2)-tilesize );
        useconv = 1;
        if( useconv == 0 )
            %Compute the distances from the template to target for all i,j
            for a = 1:size(distances,1)
                v1 = Y(startI:endI, startJ:endJ, 1:3);
                for b = 1:size(distances,2),
                    v2 = X(a:a+tilesize-1,b:b+tilesize-1, 1:3);
                    distances(a,b) = myssd( double((v1(:) > 0)) .* (v1(:) - v2(:)) );
                    %distances(a,b) = D;
                end
            end
        else
            %Compute the distances from the source to the left overlap region
            if( j > 1 )
                distances = ssd( X, Y(startI:endI, startJ:startJ+overlap-1, 1:3) );
                distances = distances(1:end, 1:end-tilesize+overlap);
            end
            %Compute the distance from the source to top overlap region
            if( i > 1 )
                Z = ssd( X, Y(startI:startI+overlap-1, startJ:endJ, 1:3) );
                Z = Z(1:end-tilesize+overlap, 1:end);
                if( j > 1 ) distances = distances + Z;
                else distances = Z;
                end
            end
            %If both are greater, compute the distance of the overlap
            if( i > 1 && j > 1 )
                Z = ssd( X, Y(startI:startI+overlap-1, startJ:startJ+overlap-1, 1:3) );
                Z = Z(1:end-tilesize+overlap, 1:end-tilesize+overlap);
                distances = distances - Z;
            end
        end
        if i==1 && j==1 && not(isempty(initial_patch))
            % initial guess
            distances = distances*0+1000;
            distances(initial_patch(1), initial_patch(2)) = 0;
        end
        %Find the best candidates for the match
        best = max( min(distances(:)), 0 );
        candidates = find(distances(:) <= (1+err)*best);
        idx = candidates(ceil(rand(1)*length(candidates)));
        [sub(1), sub(2)] = ind2sub(size(distances), idx);
        if use_verb
            fprintf( 'Picked tile (%d, %d) out of %d candidates.  Best error=%.4f\n', sub(1), sub(2), length(candidates), best );
        end
        %If we do the simple quilting (no cut), just copy image
        if( simple )
            Y(startI:endI, startJ:endJ, 1:3) = X(sub(1):sub(1)+tilesize-1, sub(2):sub(2)+tilesize-1, 1:3);
        else
            %Initialize the mask to all ones
            M = ones(tilesize, tilesize);
            %We have a left overlap
            if( j > 1 )
                %Compute the SSD in the border region
                E = ( X(sub(1):sub(1)+tilesize-1, sub(2):sub(2)+overlap-1) - Y(startI:endI, startJ:startJ+overlap-1) ).^2;
                %Compute the mincut array
                C = mincut(E, 0);
                %Compute the mask and write to the destination
                M(1:end, 1:overlap) = double(C >= 0);
            end
            %We have a top overlap
            if( i > 1 )
                %Compute the SSD in the border region
                E = ( X(sub(1):sub(1)+overlap-1, sub(2):sub(2)+tilesize-1) - Y(startI:startI+overlap-1, startJ:endJ) ).^2;
                %Compute the mincut array
                C = mincut(E, 1);
                %Compute the mask and write to the destination
                M(1:overlap, 1:end) = M(1:overlap, 1:end) .* double(C >= 0);
            end
            if( i == 1 && j == 1 )
                Y(startI:endI, startJ:endJ, 1:3) = X(sub(1):sub(1)+tilesize-1, sub(2):sub(2)+tilesize-1, 1:3);
            else
                %Write to the destination using the mask
                Y(startI:endI, startJ:endJ, :) = filtered_write(Y(startI:endI, startJ:endJ, :), ...
                    X(sub(1):sub(1)+tilesize-1, sub(2):sub(2)+tilesize-1, :), M);
            end
        end
        if use_draw
            imageplot(Y);
            drawnow;
        end
    end
end

if use_draw
    imageplot(Y);
    drawnow;
end

function y = myssd( x )
y = sum( x.^2 );

function A = filtered_write(A, B, M)
for i = 1:3,
    A(:, :, i) = A(:,:,i) .* (M == 0) + B(:,:,i) .* (M == 1);
end



%function Y = mincut(X,dir)
%Computes the minimum cut from one edge of the matrix to the other
%
%Inputs:
%   X:  Evalutations of 2d function to be cut along local minima
%   dir: 0 = vertical cut, 1 = horizontal cut
%
%Outputs:
%   C: Matrix containing entries indicating side
%       -1: left (or top) side of cut
%        0: along the cut
%        1: right (or bottom) side of cut

function C = mincut(X,dir)

if( nargin > 1 && dir == 1 )
    X = X';
end

%Allocate the current cost array, and set first row to first row of X
E = zeros(size(X));
E(1:end,:) = X(1:end,:);

%Starting with the second array, compute the path costs until the end
for i=2:size(E,1),
    E(i,1) = X(i,1) + min( E(i-1,1), E(i-1,2) );
    for j=2:size(E,2)-1,
        E(i,j) = X(i,j) + min( [E(i-1,j-1), E(i-1,j), E(i-1,j+1)] );
    end
    E(i,end) = X(i,end) + min( E(i-1,end-1), E(i-1,end) );

end

%Backtrace to find the cut
C = zeros(size(X));

[cost, idx] = min(E(end, 1:end));
C(i, 1:idx-1) = -1;
C(i, idx) = 0;
C(i, idx+1:end) = +1;

for i=size(E,1)-1:-1:1,
    for j=1:size(E,2),

        if( idx > 1 && E(i,idx-1) == min(E(i,idx-1:min(idx+1,size(E,2))) ) )
            idx = idx-1;
        elseif( idx < size(E,2) && E(i,idx+1) == min(E(i,max(idx-1,1):idx+1)) )
            idx = idx+1;
        end


        C(i, 1:idx-1) = -1;
        C(i, idx) = 0;
        C(i, idx+1:end) = +1;

    end
end

if( nargin > 1 && dir == 1 )
    %E = E';
    C = C';
end



%function Z = ssd(X, Y)
%Computes the sum of squared distances between X and Y for each possible
% overlap of Y on X.  Y is thus smaller than X
%
%Inputs:
%   X - larger image
%   Y - smaller image
%
%Outputs:
%   Each pixel of Z contains the ssd for Y overlaid on X at that pixel

function Z = ssd(X, Y)

K = ones(size(Y,1), size(Y,2));

for k=1:size(X,3),
    A = X(:,:,k);
    B = Y(:,:,k);

    a2 = filter2(K, A.^2, 'valid');
    b2 = sum(sum(B.^2));
    ab = filter2(B, A, 'valid').*2;

    if( k == 1 )
        Z = ((a2 - ab) + b2);
    else
        Z = Z + ((a2 - ab) + b2);
    end
end


