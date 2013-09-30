function [y,nbr_bits] = perform_jbig_coding(x)

% perform_jbig_coding - perform binary image coding
%
%  [y,nbr_bits] = perform_jbig_coding(x);
%
%	It requires pbmtojbg and jbgtopbm executable.
%
%	Copyright (c) 2006 Gabriel PeyrŽ

name_pbm = 'tmp.pbm';
name_jbg = 'tmp.jbg';
if size(x,1)>1 && size(x,2)>1
    % forward transform
    % save as pbm
    imwrite(rescale(x), name_pbm, 'pbm');
    % convert to jgib
    !pbmtojbg tmp.pbm tmp.jbg
    % read jbig file
    fid = fopen(name_jbg);
    if fid<0
        error('Unable to open Jbig file.');
    end
    [y,cnt] = fread(fid, Inf);
    fclose(fid);
    nbr_bits = length(y)*8;
    % remove tmp files
    !del tmp.jbg
    !del tmp.pbm
else
    % backward transform
    fid = fopen(name_jbg, 'wb');
    if fid<0
        error('Unable to open Jbig file.');
    end
    fwrite(fid, x);
    fclose(fid);
    % convert to pbm
    !jbgtopbm tmp.jbg tmp.pbm
    % read pbm
    y = imread(name_pbm);
    % remove tmp files
    !del tmp.jbg
    !del tmp.pbm
    nbr_bits = -1;
end