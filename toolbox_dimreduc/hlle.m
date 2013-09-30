function [Y, mse] = HLLE(X,k,d)

%HLLE Runs the standard Hessian LLE implementation of Hessian Eigenmaps
%   
%   X is the high-dimensional data to be processed
%   k is the number of nearest neighbor points to be used 
%   if k is a scalar, same size used at all points. if k is a vector of
%   length N, neighborhood N(i) will be assigned k(i) nearest neighbors
%   
%   d is the number of dimensions to embed X in
%   
%   Y is the output embedded data
%   mse is the sum (at each neighborhood used) of the eigenvalues(d+2:end)
%   of the local coordinate representation. used for adaptive neighborhood
%   restructuring
%   Example:
%   N=1000; k=12; d=2;
%   tt = (3*pi/2)*(1+2*rand(1,N));  height = 21*rand(1,N);
%   X = [tt.*cos(tt); height; tt.*sin(tt)];
%   [Y, mse] = HLLE(X,k,d);

%   C. Grimes and D. Donoho, March 2003
%   Last Revision: 
%   Version 1.0

%get data size
N = size(X,2);
%check for constant neighborhood size
if max(size(k)) ==1
    kvec = repmat(k,N,1);
elseif max(size(k)) == N
    kvec=k;
else
    error('Neighborhood Vector Size does not match data');
end;

disp(['-->Running HLLE for ', num2str(N), ' points']);
disp('-->Computing HLLE neighbors');

%Compute Nearest neighbors
D1 = L2_distance(X,X,1);

dim = size(X,1);
nind = repmat(0, size(D1,1), size(D1,2));
%extra term count for quadratic form
dp = d*(d+1)/2;
W = repmat(0,dp*N,N);

if(mean(k)>d) 
  disp('[note: k>d; regularization will be used]'); 
  tol=1e-3; % regularlizer in case constrained fits are ill conditioned
else
  tol=0;
end;

for i=1:N
    tmp = D1(:,i);
    [ts, or] = sort(tmp);
%take k nearest neighbors
    nind(or(2:kvec(i)+1),i) = 1;
    thisx = X(:,or(2:kvec(i)+1));
    %center using the mean 
    thisx = thisx - repmat(mean(thisx')',1,kvec(i));

    %compute local coordinates
    [U,D,Vpr] = svd(thisx);
    V = Vpr(:,1:d);
    
    %Neighborhood diagnostics
    vals = diag(D);
    mse(i) = sum(vals(d+1:end));
    
    
%build Hessian estimator
    clear Yi; clear Pii;
    ct = 0;
    for mm=1:d
        startp = V(:,mm);
        for nn=1:length(mm:d)
            indles = mm:d;
            Yi(:,ct+nn) = startp.*(V(:,indles(nn)));
        end;
        ct = ct+length(mm:d);
    end;
    Yi = [repmat(1,kvec(i),1), V, Yi];
%orthogonalize linear and quadratic forms
    [Yt, Orig] = mgs(Yi);
    Pii = Yt(:,d+2:end)';
 %double check weights sum to 1
    for j=1:dp
        if sum(Pii(j,:)) >0.0001
            tpp = Pii(j,:)./sum(Pii(j,:)); 
        else
            tpp = Pii(j,:);
        end;
        %fill weight matrix
       W((i-1)*dp+j, or(2:kvec(i)+1)) = tpp;
    end;
end;

%%%%%%%%%%%%%%%%%%%%Compute eigenanalysis of W
disp('-->Computing HLLE embedding');

G=W'*W;
G = sparse(G);

options.disp = 0; 
options.isreal = 1; 
options.issym = 1;

%tol=1e-3; %sometimes useful for pathological fits
tol=0;
% [Yo,eigenvals] = eigs(G,d+1,tol,options); % old code
[Yo,eigenvals] = eigs(G,d+1,'SM',options);
Y = Yo(:,1:d)'*sqrt(N); % bottom evect is [1,1,1,1...] with eval 0


disp('-->Orienting Coordinates');
%compute final coordinate alignment
R = Y'*Y;
R2 = R^(-1/2);
Y = Y*R2;


function d = L2_distance(a,b,df)
% L2_DISTANCE - computes Euclidean distance matrix
%
% E = L2_distance(A,B)
%
%    A - (DxM) matrix 
%    B - (DxN) matrix
%    df = 1, force diagonals to be zero; 0 (default), do not force
% 
% Returns:
%    E - (MxN) Euclidean distances between vectors in A and B
%
%
% Description : 
%    This fully vectorized (VERY FAST!) m-file computes the 
%    Euclidean distance between two vectors by:
%
%                 ||A-B|| = sqrt ( ||A||^2 + ||B||^2 - 2*A.B )
%
% Example : 
%    A = rand(400,100); B = rand(400,200);
%    d = distance(A,B);

% Author   : Roland Bunschoten
%            University of Amsterdam
%            Intelligent Autonomous Systems (IAS) group
%            Kruislaan 403  1098 SJ Amsterdam
%            tel.(+31)20-5257524
%            bunschot@wins.uva.nl
% Last Rev : Wed Oct 20 08:58:08 MET DST 1999
% Tested   : PC Matlab v5.2 and Solaris Matlab v5.3

% Copyright notice: You are free to modify, extend and distribute 
%    this code granted that the author of the original code is 
%    mentioned as the original author of the code.

% Fixed by JBT (3/18/00) to work for 1-dimensional vectors
% and to warn for imaginary numbers.  Also ensures that 
% output is all real, and allows the option of forcing diagonals to
% be zero.  

if (nargin < 2)
   error('Not enough input arguments');
end

if (nargin < 3)
   df = 0;    % by default, do not force 0 on the diagonal
end

if (size(a,1) ~= size(b,1))
   error('A and B should be of same dimensionality');
end

if ~(isreal(a)*isreal(b))
   disp('Warning: running distance.m with imaginary numbers.  Results may be off.'); 
end

if (size(a,1) == 1)
  a = [a; zeros(1,size(a,2))]; 
  b = [b; zeros(1,size(b,2))]; 
end

aa=sum(a.*a); bb=sum(b.*b); ab=a'*b; 
d = sqrt(repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab);

% make sure result is all real
d = real(d); 

% force 0 on the diagonal? 
if (df==1)
  d = d.*(1-eye(size(d)));
end

function [Q, R] = mgs(A);
%
%   modified Gram-Schmidt. this is a more stable way to compute a
%   qr factorization
%
   [m, n] = size(A);
%
%   we assume that m>= n.
%   
   V = A;
   R = zeros(n,n);
   for i=1:n
      R(i,i) = norm(V(:,i));
      V(:,i) = V(:,i)/R(i,i);
      if (i < n)
         for j = i+1:n
            R(i,j) = V(:,i)' * V(:,j);
            V(:,j) = V(:,j) - R(i,j) * V(:,i);
         end
      end
   end
   Q = V;
