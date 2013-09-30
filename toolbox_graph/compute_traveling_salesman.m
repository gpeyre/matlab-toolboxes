function I = compute_traveling_salesman(W, options)

% compute_traveling_salesman - wrapper for the linkern programm solving TSP
%
%   I = wrapper_linkern(W, options);
%
%   W is the matrix of weights (should be intergers).
%   I is the ordering solution of the traveling salesman problem.
%
%   linkern program can be downloaded from 
%   http://www.tsp.gatech.edu/concorde/downloads/downloads.htm
%   Linkern is an implementation of the Chained-Lin-Kernighan heuristic for
%   the TSP. It solves *approximatly* the TSP in fast time.
%
%   Copyright (c) 2006 Gabriel Peyré

options.null = 0;
if isfield(options, 'nbr_runs')
    nbr_runs = options.nbr_runs;
else
    nbr_runs = 1;
end
if isfield(options, 'kick')
    kick = options.kick;
else
    kick = 3;
end

tmpfile = 'tmp.tsp';

fid = fopen(tmpfile,'W');
if fid<0
    error('Unable to open tmp file.');
    return;
end

n = size(W,1);

if max( abs(round(W(:))-W(:)) )>eps
    % convert to integers
    W = round(W*1000);
end

% write data to file
fprintf(fid, ['NAME : tmp\nTYPE : TSP\n' ...
    'COMMENT : tmp file from wrapper_linkern.\n' ...
    'DIMENSION : ' num2str(n) '\n' ...
    'EDGE_WEIGHT_TYPE : EXPLICIT\n' ...
    'EDGE_WEIGHT_FORMAT : FULL_MATRIX\n' ... 
    'EDGE_WEIGHT_SECTION\n'] );

fprintf(fid, '%d ', W(:)');
fprintf(fid, '\nEOF\n');
fclose(fid);

% call linker
!linkern -o tmp.res tmp.tsp > tmp.out
% system(['linkern -o tmp.res -r ' num2str(nbr_runs) ' -K ' num2str(kick) ' tmp.tsp > tmp.out']);

% retrieve result
fid = fopen('tmp.res','r');
if fid<0
    error('Unable to open tmp file.');
    return;
end
[A,cnt] = fscanf(fid,'%d %d', 2);
[A,cnt] = fscanf(fid,'%d %d %d');
if cnt~=3*n
    warning('Problem in reading results.');
end
A = reshape(A, 3, cnt/3);
I = A(1,:) + 1;
fclose(fid);

delete('tmp.res');
delete(tmpfile);
delete('tmp.out');