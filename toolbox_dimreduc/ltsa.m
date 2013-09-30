% ltsa - local tangent planes alignement
%
%   [T,NI] = ltsa(X,d,K);
%
%   X is the (d,n) n data points in R^d.
%   d is the output dimensionnality.
%   K is the number of nearest neighbors.
%
% Written by Zhenyue Zhang & Hongyuan Zha, 2004.
% Reference: http://epubs.siam.org/sam-bin/dbq/article/41915

function [T,NI] = ltsa(X,d,K,NI)

[m,N] = size(X);  % m is the dimensionality of the input sample points.
% Step 0:  Neighborhood Index
if nargin<4
    if length(K)==1
        K = repmat(K,[1,N]);
    end;
    NI = cell(1,N);
    if m>N
        if 0
            a = sum(X.*X);
            dist2 = sqrt(repmat(a',[1 N]) + repmat(a,[N 1]) - 2*(X'*X));
            for i=1:N
                % Determine ki nearest neighbors of x_j
                [dist_sort,J] = sort(dist2(:,i));
                Ii = J(1:K(i));
                NI{i} = Ii;
            end;
        else
            % use fast code
            options.exlude_self = 1;
            [D1,nn_list] = compute_nn_distance(X,K(1), options);
            for i=1:N
                NI{i} = nn_list(i,:);
            end
        end
    else
        if 0
            for i=1:N
                % Determine ki nearest neighbors of x_j
                x = X(:,i); ki = K(i);
                dist2 = sum((X-repmat(x,[1 N])).^2,1);
                [dist_sort,J] = sort(dist2);
                Ii = J(1:ki);
                NI{i} = Ii;
            end;
        else
            % use fast code
            options.exlude_self = 1;
            [D1,nn_list] = compute_nn_distance(X,K(1), options);
            for i=1:N
                NI{i} = nn_list(i,:);
            end
        end
    end;
else
    K = zeros(1,N);
    for i=1:N
        K(i) = length(NI{i});
    end;
end;
% Step 1:  local information
BI = {};
Thera = {};
for i=1:N
    % Compute the d largest right singular eigenvectors of the centered matrix
    Ii = NI{i}; ki = K(i);
    Xi = X(:,Ii)-repmat(mean(X(:,Ii),2),[1,ki]);
    W = Xi'*Xi; W = (W+W')/2;
    [Vi,Si] = schur(W);
    [s,Ji] = sort(-diag(Si));
    Vi = Vi(:,Ji(1:d));
    % construct Gi
    Gi = [repmat(1/sqrt(ki),[ki,1]) Vi];
    % compute the local orthogonal projection Bi = I-Gi*Gi'
    % that has the null space span([e,Theta_i^T]).
    BI{i} = eye(ki)-Gi*Gi';
end;
B = speye(N);
for i=1:N
    Ii = NI{i};
    B(Ii,Ii) = B(Ii,Ii)+BI{i};
    B(i,i) = B(i,i)-1;
end;
B = (B+B')/2;
options.disp = 0;
options.isreal = 1;
options.issym = 1;
[U,D] = eigs(B,d+2,0,options);
lambda = diag(D);
[lambda_s,J] = sort(abs(lambda));
U = U(:,J); lambda = lambda(J);
T = U(:,2:d+1)';

