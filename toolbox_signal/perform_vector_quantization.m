function [d,e] = perform_vector_quantization(a,b,c)

% perform_vector_quantization - perform the quantization of a matrix
%
%   Forward (quantization)
% [C,IDX] = perform_vector_quantization(X,options);
%   Backward (de-quantization)
% X = perform_vector_quantization(C,IDX,options);
%
%   X is a 1D/2D/3D matrix.
%   options.block_size is the size of the block used for creating the quantized vectors.
%   options.nbr_iterations is the number of iteration of the Lloyd smoothing for the computation
%       of the dictionary.
%   options.nbr_vectors is the number of vector in the dictionnary.
%   options.nbr_samples is the number of samples used to compute the
%       dictionnary (set it < to the number of vectors to speed up
%       computations).
%
%   Copyright (c) 2005 Gabriel PeyrŽ

options.null = 0;
if nargin==1
    dir = 1;
    X = a;
elseif nargin==2 && length(b)<=1
    dir = 1;
    X = a;
    options = b;
elseif nargin==2
    dir = -1;
    C = a;
    IDX = b;
else
    dir = -1;
    C = a; IDX = b; options = c;
end
    
% size of the blocks
if isfield(options, 'block_size')
    Z = options.block_size;
else
    Z = [3 3 2];
end

% dimension of the vectors
d = prod(Z);
% to be sure to have block of length 3
Z(end+1:3) = 1;

if dir==-1
    % bwd transform
    nbr_vectors = length(IDX);  % number of vectors
    dim = size(C,1);            % dimension of vectors
    if dim~=prod(Z)
        error('Size of vectors does not match size of blocks');
    end
    % look up
    M = C(:,IDX(:));
    % reshaping
    s = size(IDX); s(end+1:3) = 1;
    A = Z .* s;
    X = zeros( A );
    m = 0;
    for i=1:size(IDX,1)
        for j=1:size(IDX,2)
            for k=1:size(IDX,3)
                m = m+1;
                selx = (i-1)*Z(1)+1:i*Z(1);
                sely = (j-1)*Z(2)+1:j*Z(2);
                selz = (k-1)*Z(3)+1:k*Z(3);
                X(selx,sely,selz) = reshape( M(:,m), Z(1), Z(2), Z(3) );
            end
        end
    end
    d = X;
    return;
end

if isfield(options, 'verb')
    verb = options.verb;
else
    verb = 1;
end

% nbr of iteration for Lloyd clustering
if isfield(options, 'nbr_iterations')
    nbr_iterations = options.nbr_iterations;
else
    nbr_iterations = 10;
end



% extend by zeros to get integer size
s = size(X); s(end+1:3) = 1;
a = ceil( s./Z );
s1 = a.*Z;
X1 = zeros(s1);
X1(1:s(1),1:s(2),1:s(3)) = X;
X = X1; clear X1;

% total number of points
p = prod(a);
V = zeros(d,p);

% number of vectors
if isfield(options, 'nbr_vectors')
    nbr_vectors = options.nbr_vectors;
else
    nbr_vectors = min(128, p/4);
end

m = 0;
for i = 1:a(1)
    for j = 1:a(2)
        for k = 1:a(3)
            m = m+1;
            selx = (i-1)*Z(1)+1:i*Z(1);
            sely = (j-1)*Z(2)+1:j*Z(2);
            selz = (k-1)*Z(3)+1:k*Z(3);
            v = X(selx,sely,selz);
            V(:,m) = v(:);
        end
    end
end

% use only a small subset to speed up
if isfield(options, 'nbr_samples')
    nbr_samples = options.nbr_samples;
else
    nbr_samples = min(1024*2, p);
end

nbr_samples = min(nbr_samples,p);

% select at random
s = rand(p,1);
[tmp,s] = sort(s);
s = sort( s(1:nbr_samples) );

% perform lloyd iteration
if verb
    fprintf('Performing Lloyd iteration ... ');
end
options.nb_iter = nbr_iterations;
[IDX,C] = perform_kmeans(V(:,s)',nbr_vectors,options);
% [IDX,C] = kmeans(V(:,s)',nbr_vectors, 'Maxiter', nbr_iterations, 'EmptyAction', 'singleton'); 
% C = C';
if verb
    fprintf('done.\n');
end

% perform clustering
if verb
    fprintf('Performing clustering ... ');
end
D = compute_distance_to_points(V,C); 
[tmp,IDX] = min(D); IDX = IDX(:);
IDX = reshape(IDX,a(1),a(2),a(3));
if verb
    fprintf('done.\n');
end

% assign result
d = C; e = IDX;