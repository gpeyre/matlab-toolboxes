function write_lum(M,name)

% write_lum - write a file with .lum extension
%
%   write_lum(M,name);
%
%   The lum format is usualy a 12 bit file format, written as
%   a 32 bit float.
%
%   Copyright (c) 2006 Gabriel Peyré

do_rounding = 0;

if do_rounding
    M = round(M);
end

fid = fopen(name, 'wb');
if fid<0
    error(['Impossible to open file ' name ' for writing.']);
end 
[n,p] = size(M);
% first write size
v = zeros(n,1); v(1) = n; v(2) = p;
fwrite(fid, v(:), 'float32');
% then write content
fwrite(fid, M(:), 'float32');
fclose(fid);
