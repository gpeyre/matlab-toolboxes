function D = perform_dijkstra_propagation_old( G , S )

% dijkstra - Find shortest paths in graphs
%
% 	D = dijkstra_fast( G , S );
%
%   use the full or sparse matrix G in which 
% 	an entry (i,j) represents the arc length between nodes i and j in a 
%	graph. In a full matrix, the value INF represents the absence of an arc; 
%	in a sparse matrix, no entry at (i,j) naturally represents no arc.
%		
%	S is the one-dimensional matrix of source nodes for which the shortest
%	to ALL other nodes in the graphs will be calculated. The output matrices
%	D and P contain the shortest distances and predecessor indices respectively.
%	An infinite distance is represented by INF. The predecessor indices contain
%	the node indices of the node along the shortest path before the destination
%	is reached. These indices are useful to construct the shortest path with the
%	function pred2path (by Michael G. Kay).
%
%	This function was implemented in C++. The source code can be compiled to a
%	Matlab compatible mex file by the command "mex -O dijkstra.cpp" at the Matlab
%	prompt. In this package, we provide a compiled .dll version that is 
%       compatible all Windows based machines.  If you are not working on a 
%       Windows platform, delete the .dll version provided and recompile from
%       the .cpp source file.  If you do not have the Matlab compiler or a Windows
%       platform, delete the .dll version and dijkstra will then call the
%	Matlab function dijk.m (by Michael G. Kay).  Note that this Matlab
%	code is several orders of magnitude slower than the C based mex file.
%
% code taken from  :
%   Mark Steyvers, Stanford University, 12/19/00

N = size( G , 1 );
D = dijk( G , S , 1:N );





function D = dijk(A,s,t)
% dijk - shortest paths from nodes 's' to nodes 't' using Dijkstra algorithm.
%
%   D = dijk(A,s,t);
%
%     A = n x n node-node weighted adjacency matrix of arc lengths
%         (Note: A(i,j) = 0   => Arc (i,j) does not exist;
%                A(i,j) = NaN => Arc (i,j) exists with 0 weight)
%     s = FROM node indices
%       = [] (default), paths from all nodes
%     t = TO node indices
%       = [] (default), paths to all nodes
%     D = |s| x |t| matrix of shortest path distances from 's' to 't'
%       = [D(i,j)], where D(i,j) = distance from node 'i' to node 'j' 
%
%	(If A is a triangular matrix, then computationally intensive node
%   selection step not needed since graph is acyclic (triangularity is a 
%   sufficient, but not a necessary, condition for a graph to be acyclic)
%   and A can have non-negative elements)
%
%	(If |s| >> |t|, then DIJK is faster if DIJK(A',t,s) used, where D is now
%   transposed and P now represents successor indices)
%
%  (Based on Fig. 4.6 in Ahuja, Magnanti, and Orlin, Network Flows,
%   Prentice-Hall, 1993, p. 109.)

% Copyright (c) 1998-2000 by Michael G. Kay
% Matlog Version 1.3 29-Aug-2000
% 
%  Modified by JBT, Dec 2000, to delete paths


% Input Error Checking ******************************************************
error(nargchk(1,3,nargin));

[n,cA] = size(A);

if nargin < 2 | isempty(s), s = (1:n)'; else s = s(:); end
if nargin < 3 | isempty(t), t = (1:n)'; else t = t(:); end

if ~any(any(tril(A) ~= 0))			% A is upper triangular
   isAcyclic = 1;
elseif ~any(any(triu(A) ~= 0))	% A is lower triangular
   isAcyclic = 2;
else										% Graph may not be acyclic
   isAcyclic = 0;
end

if n ~= cA
   error('A must be a square matrix');
elseif ~isAcyclic & any(any(A < 0))
   error('A must be non-negative');
elseif any(s < 1 | s > n)
   error(['''s'' must be an integer between 1 and ',num2str(n)]);
elseif any(t < 1 | t > n)
   error(['''t'' must be an integer between 1 and ',num2str(n)]);
end
% End (Input Error Checking) ************************************************

A = A';		% Use transpose to speed-up FIND for sparse A

D = zeros(length(s),length(t));
P = zeros(length(s),n);

for i = 1:length(s)
   j = s(i);
   
   Di = Inf*ones(n,1); Di(j) = 0;
   
   isLab = logical(zeros(length(t),1));
   if isAcyclic ==  1
      nLab = j - 1;
   elseif isAcyclic == 2
      nLab = n - j;
   else
      nLab = 0;
      UnLab = 1:n;
      isUnLab = logical(ones(n,1));
   end
   
   while nLab < n & ~all(isLab)
      if isAcyclic
         Dj = Di(j);
      else	% Node selection
         [Dj,jj] = min(Di(isUnLab));
         j = UnLab(jj);
         UnLab(jj) = [];
         isUnLab(j) = 0;
      end
      
      nLab = nLab + 1;
      if length(t) < n, isLab = isLab | (j == t); end
      
      [jA,kA,Aj] = find(A(:,j));
      Aj(isnan(Aj)) = 0;
            
      if isempty(Aj), Dk = Inf; else Dk = Dj + Aj; end
      
      P(i,jA(Dk < Di(jA))) = j;
      Di(jA) = min(Di(jA),Dk);
      
      if isAcyclic == 1			% Increment node index for upper triangular A
         j = j + 1;
      elseif isAcyclic == 2	% Decrement node index for lower triangular A
         j = j - 1;
      end
      
      %disp( num2str( nLab ));
   end
   D(i,:) = Di(t)';
end

