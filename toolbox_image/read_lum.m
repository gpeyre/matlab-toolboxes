function M = read_lum(name)

% read_lum - read a file with .lum extension
%
%   M = read_lum(name);
%
%   The lum format is usualy a 12 bit file format, written as
%   a 32 bit float.
%
%   Copyright (c) 2006 Gabriel Peyré

fid = fopen(name, 'rb');
if fid<0
    error(['File ' name ' does not exists.']);
end 
A = fread(fid, Inf, 'float32');
fclose(fid);

% test if the size of the file is written
n = A(1); p = A(2);
if length(A)==(p+1)*n
    A(1:n)=[];
    M = reshape(A, n,p);
else
    n = sqrt(length(A)-3); p = n;
    if length(A)~=n^2+3
        error('Problem during conversion (e.g. image is not square.)');
    end
    M = reshape(A(4:end), n,p);
end