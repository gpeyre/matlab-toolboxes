% perform_dijkstra_propagation - shortest distance on graph
%
%   This is impemented in a mex file.
%
%   [D,S] = perform_dijkstra_propagation(W,start_verts,end_verts,nb_iter_max,H);
%    D is the distance to starting points.
%    S is the state : dead=-1, open=0, far=1.
%    
%    Copyright (c) 2005 Gabriel Peyré