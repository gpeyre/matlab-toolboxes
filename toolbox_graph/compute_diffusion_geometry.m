function V = compute_diffusion_geometry(T,precision)

% compute_diffusion_geometry - compute the diffusion geometry basis functions
%
%   V = compute_diffusion_geometry(T,precision);
%
%   T is a diffusion operator associated.
%   V is the wavelet orthogonal basis associated to this operator.
%
%   See 
%       R.R.Coifman, S. Lafon, A.B. Lee, M. Maggioni, B. Nadler, F. Warner and S.W. Zucker, 
%       “Geometric diffusions as a tool for harmonic analysis and structure definition of data”, 
%       Proceedings of the National Academy of Sciences (May 2005)
%
%   This function uses the code from Mauro Maggioni
%       <http://www.math.yale.edu/%7Emmm82/diffusionwavelets.html>
%   that should be in your path
%
%   Copyright (c) 2005 Gabriel Peyré

if nargin<2
    precision = 1e-10;
end

n = size(T,1);


% add the toolbox of Maggioni in your path
if exist('diffusion_wavelets')>0
    path(path, 'diffusion_wavelets/');
end

% use default settings to compute the wavelet basis
Levels = 10;
Type = 'Full';
Type = 'Partial';
TreeOptions = [];
Tree = DiffusionWaveletTree( Type, TreeOptions, speye(n), T, Levels, precision);

V = zeros(n);

% retrieve the vectors one by one.
k = 0;
for i = 2:length(Tree)
    s = 2;
    nbr = size( Tree{i}(s).Basis, 2 );
    for num = 1:nbr
        k = k+1;
        F = GetBasisFcn(Tree, i, s, num);
        V(:,k) = F;
    end
end

% retrieve low frequencies
i = length(Tree); s = 1;
nbr = size( Tree{i}(s).Basis, 2 );
for num = 1:nbr
    k = k+1;
    F = GetBasisFcn(Tree, i, s, num);
    V(:,k) = F;
end

% reverse so that low frequencies come first
V = V(:,end:-1:1);