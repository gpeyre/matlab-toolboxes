function T = perform_voting(M,sigma,c,n,p,thresh)

% perform_voting - perform the voting process.
%
%   T = perform_voting(M,sigma,c,n,p,thresh);
%
%   'M' is a tensor field created by the user.
%
%   - The oriented features of the image must be stored
%     in the *highest* eigenvectors of the tensors.
%   - The confidence in these featres is reflected by
%     ratio of the eigenvalues 
%     (l1/l2==1 imply no confidence at all, l2=0 imply total confidence).
%
%   Optional:
%   'sigma' and 'c' control the stick field (see 'compute_stick_tf').
%   'n' is the size of the kernel (typically smaller than the size of the image).
%   'p' is the precision for the computation of the ball tensor field.
%   'thresh' is a threshold (perform voting only if l1-l2>thresh where (l1,l2) are the eigenvalues).
%
%   See also perform_voting_ball.
%
%   Copyright (c) 2004 Gabriel Peyré

verb = 1;

N = size(M,1);
if nargin<2
    sigma = 0.2;
end
if nargin<3
    c = 0.1*sigma;
end
if nargin<4
    n = N/2;
end
if nargin<5
    p = 32;
end
n = floor(n/2)*2+1;
nn = (n-1)/2;

% perform eigen-decomposition
[e1,e2,l1,l2] = perform_tensor_decomp(M);
    
% perform ball voting
T = perform_voting_ball(l2,sigma,c);

% perform a threshold to use only high coef
if nargin<6
    thresh = mean(mean( (l1-l2) ))/10;
end
I = find( l1-l2>thresh );

% perform stick voting
if verb
    h = waitbar(0,'Performing voting...');
end
a = 0;
for s = I';
    a = a+1;
    if verb
        waitbar(a/length(I),h);
    end
    [i,j] = ind2sub(size(l1),s);
    % the direction is e1 with intensity l1-l2
    v = e1(i,j,:);    
    Fk = compute_stick_tf(v(:),N,n,sigma,c);    
    
    if 0
        [e1,e2,l1,l2] = perform_tensor_decomp(Fk);
        clf;
        plot_vf(e1,l1);
    end
    
    Gk = l1(i,j)-l2(i,j);    % the current intensity of the vote
    % select the correct range
    i1 = max(1-(i-nn),0);
    i2 = max(i+nn-N,0);
    j1 = max(1-(j-nn),0);
    j2 = max(j+nn-N,0);
    T( i-nn+i1:i+nn-i2, j-nn+j1:j+nn-j2,:,: ) = T( i-nn+i1:i+nn-i2, j-nn+j1:j+nn-j2,:,: ) + Gk*Fk(1+i1:end-i2, 1+j1:end-j2,:,:);
end
if verb
    close(h);
end