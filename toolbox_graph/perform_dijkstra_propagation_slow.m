function [D,S] = perform_dijkstra_propagation_slow(W,start_verts,end_verts,nb_iter_max,H,verbose)

% perform_dijkstra_propagation_slow - Matlab implementation of Dijkstra algorithm.
%
%   [D,S] = perform_dijkstra_propagation_slow(W,start_verts,end_verts,nb_iter_max,H);
% 
%	D is the distance to starting points.
%	S is the state : dead=-1, open=0, far=1.
%
%   Copyright (c) 2005 Gabriel Peyré


%   'data' is a structure containing all data for the dijkstra algorithm:
%       - 'data.A': action (distance to starting points)
%       - 'data.O': open list
%       - 'data.C': close list
%       - 'data.S': state, either 'O' or 'C'
%       - 'data.F': Father
%       - 'data.P': Origin seed point
%       - 'data.H': Heuristic (when the dijkstra is used as A*)


if nargin<3
    end_verts = [];
end
if nargin<4
    nb_iter_max = Inf;
end
if nargin<6
	verbose = 1;
end

n = size(W,1);
if sum(start_verts>n)>0 || sum(end_verts>n)>0
    error('Out of bound start/end vertices.');
end


nb_iter_max = min(nb_iter_max, size(W,1));

% initialize the data
data = dijkstra_init(W, start_verts);

str = 'Performing Dijkstra algorithm.';
if verbose
    h = waitbar(0,str);
end

i = 0; 
while i<nb_iter_max
    i = i+1;
    data = dijkstra_step(data);
    if verbose
        waitbar(i/nb_iter_max, h, sprintf('Performing Dijkstra algorithm, iteration %d.', i) );
     end
    % check if we have reach one of the end points
    for j=end_verts
        if ~isempty( find(data.C==j) )
            if verbose
                close(h);
            end
            S = data.S;
            D = data.A;
            I = find(S==0); S(I) = 1;
            I = find(S=='O'); S(I) = 0;
            I = find(S=='C'); S(I) = -1;
            return;
        end
    end
end

if verbose
    close(h);
end

S = data.S;
D = data.A;
I = find(S=='O'); S(I) = 0;
I = find(S==0); S(I) = 1;
I = find(S=='D'); S(I) = -1;

function data = dijkstra_init(W, start_verts, heuristic)

% dijkstra_init - initialisation of dijkstra algorithm
%
%   data = dijkstra_init(W, start_verts [,heuristic]);
%
%   'heuristic' is a structure that should contains :
%       - one field 'heuristic.func' which should be a function.
%           This function take as input the number of a vertex
%           and should return a heuristical measure of the distance
%           between this point and the target.
%       - one field 'heuristic.weight' in [0,1] which measure
%           how much this heuristic should be taken into acount.
%           0 is classical Dijkstra, and 1 is full A* algorithm.
%       - you can add other fields (use data) for your function.
%
%   Copyright (c) 2004 Gabriel Peyré

n = size(W,1);

if nargin<3
    heuristic.func  = 0;
    heuristic.weight  = 0;
end

data.heuristic = heuristic;

data.A = zeros(n,1) + Inf; % action 
data.A(start_verts) = 0;

data.O = start_verts;

data.C = [];

data.F = zeros(n,1) - 1;
data.P = zeros(n,1) - 1;
data.P(start_verts) = start_verts;
data.H = zeros(n,1);
data.S = zeros(n,1);
data.S(start_verts) = 'O';

data.adj_list = adjmatrix2list(W);

data.W = W;


function data1 = dijkstra_step(data)

% dijkstra_step - perform one step in the dijkstra algorithm
%
%   [O1,C1] = dijkstra_step(O,C,W,adj_list);
%
%   Copyright (c) 2004 Gabriel Peyré

A = data.A; % action 
O = data.O; % open list
C = data.C; % close list
S = data.S; % state, either 'O' or 'C'
F = data.F; % Father
P = data.P; % Origin seed point
H = data.H; % Heuristic
adj_list = data.adj_list;   % adjacency list
W = data.W; % weight matrix
heuristic = data.heuristic;

if isempty(O)
    data1 = data;
    return;
end

[m,I] = min(A(O)+H(O));
x = O(I(1));   % selected vertex

% pop from open and add to close
O = O( find(O~=x) );
C = [C,x];
S(x) = 'C'; % now its close
% its neighbor
nei = adj_list{x}; 

for i=nei
    w = W(x,i);
    A1 = A(x) + w;    % new action from x
    switch S(i)
        case 'C'
            % check if action has change. Should not appen for dijkstra
            if A1<A(i)
                % pop from Close and add to Open  
                C = C( find(C~=i) );
                O = [O,i];
                S(i) = 'O';
                A(i) = A1;
                F(i) = x;       % new father in path
                P(i) = P(x);    % new origin
            end
        case 'O'
            % check if action has change.
            if A1<A(i)
                A(i) = A1;
                F(i) = x;   % new father in path
                P(i) = P(x);    % new origin
            end
        otherwise
            if A(i)~=Inf
                warning('Action must be initialized to Inf');  
            end    
            if heuristic.weight~=0
                % we use an heuristic
                str = sprintf( 'd = %s( i, heuristic );', heuristic.func );
                eval(str);
                H(i) = heuristic.weight*d;
            end
            % add to open
            O = [O,i];
            S(i) = 'O';
            % action must have change.
            A(i) = A1;
            F(i) = x;   % new father in path
            P(i) = P(x);    % new origin
    end
end

data1.A = A;
data1.O = O;
data1.C = C;
data1.S = S;
data1.P = P;
data1.adj_list = adj_list;
data1.W = W;
data1.F = F;
data1.heuristic = heuristic;
data1.H = H;