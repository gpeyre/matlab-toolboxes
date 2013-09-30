function [y,nbr_bits] = perform_jp2k_compression(x,options)

% perform_jp2k_compression - compress/decompress an image using JPEG2000
%
% Compression:
%   [stream,nbr_bits] = perform_jp2k_compression(M,options);
% De-compression:
%   M = perform_jp2k_compression(stream,options);
%
%   M is the image to code (2D array or 3D for colors).
%   It should be integer valued (support for up 16 bits per pixel).
%
%   This function use OpenJPEG implementation of JPEG2000.
%   http://www.tele.ucl.ac.be/PROJECTS/OPENJPEG/
%
%   You must have in your path the binary executables
%   image_to_j2k and j2k_to_image.
%
%   Copyright (c) 2006 Gabriel Peyré

if size(x,1)==1 || size(x,2)==1
    dir = -1;
else
    dir = 1;
end

if isfield(options, 'db')
    db = [' -q ' num2str(options.db)];
else
    db = [];
end
if isfield(options, 'rate')
    rate = [' -r ' num2str(options.rate)];
else
    rate = [];
end
if isfield(options, 'nbr_bits')
    if dir==1
        nbr_bits = options.nbr_bits;
        % transform into rate
        if ~isempty(rate)
            warning('Cannot use both options.db and options.rate');
        else
            bd = ceil(log2(max(x(:))));
            nbr_bits_orig = bd*prod(size(x));
            rate = round(nbr_bits_orig/nbr_bits);
            rate = [' -r ' num2str(rate)];
        end
    end
end

if dir==1
    if max(x(:))>255 || min(x(:))<0
        warning('Out of range data.');
    end
    % write the image as 16 bit pgm format
    warning off;
    xc = clamp(x/255,0,1);
    imwrite(xc, 'tmp.pgm');
    warning on;
    % perform compression
    system(['./image_to_j2k -i tmp.pgm -o tmp.j2k -I' db rate]);
    % read back from file
    fid = fopen('tmp.j2k');
    if fid<0
        error('Unable to open j2k file.');
    end
    [y,cnt] = fread(fid, Inf);
    fclose(fid);
    nbr_bits = length(y)*8;
else
    % write into file
    fid = fopen('tmp.j2k', 'wb');
    if fid<0
        error('Unable to open j2k file.');
    end
    fwrite(fid, x);
    fclose(fid);
    % decompress
    !j2k_to_image -i tmp.j2k -o tmp.pgm
    y = double( imread('tmp.pgm') );
end

% delete tmp files
warning off;
delete tmp.*
warning on;