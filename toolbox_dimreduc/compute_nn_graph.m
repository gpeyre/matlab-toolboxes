function [D,A] = compute_nn_graph(X,options)

% compute_nn_graph - compute nearest neighbor graph
%
% [D,A] = compute_nn_graph(X,options);
%
%	X is an array of points, size must be (dimension x nbr_points)
%       X(:,i) is the ith point.
%
%   D is the graph matrix : D(i,j)==Inf if point i and j are not connected,
%       D(i,j)=|X(:,i)-X(:,j)| otherwise
%
% The graph can be computed using two modes:
%   options.nn_epsilon=dist_between_NN  ==> use epsilon-based NN mode
%   options.nn_nbr=#NN                  ==> use k-based NN mode.
%   options.nn_automatic=#NN            ==> evaluate epsilon to match the
%           k-NN (tends to prune outliers.
%
%   Copyright (c) 2005 Gabriel Peyré


if isfield(options, 'verb')
    verb = options.verb;
else
    verb = 1;
end


if isfield(options, 'nn_epsilon')
    nn_epsilon = options.nn_epsilon;
    nn_mode = 'epsilon';
else
    nn_mode = 'nbr';
    if isfield(options, 'nn_nbr')
        nn_nbr = options.nn_nbr;
    else
        nn_nbr = 5;
    end
end




if strcmp(nn_mode, 'epsilon')
    if verb
        disp('- computing distance matrix.');
    end
    D = sqrt( compute_distance_matrix(X) );
    % mode by thresholding
    I = find(D>nn_epsilon);
    D(I) = Inf;
else
    % D = sqrt( compute_distance_matrix(X) );
    if verb
        disp('- computing nn graph.');
    end
    options.exlude_self = 1;
    [D1,nn_list] = compute_nn_distance(X,nn_nbr, options);
    n = size(D1,1);
    D = zeros(n)+Inf;
    % mode by selection
    for i=1:size(D,1)
        D( i,nn_list(i,:) ) = D1(i,:);
    end
end

A = ~isinf(D);

% ensure D is symmetric
D = min(D,D');