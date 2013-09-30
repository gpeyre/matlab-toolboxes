function [I_new] = jp2_degrade(Isrc, comp_ratio, bit_depth)
% Original jp2k compression/decompression function
%
% Attention: a) no rounding operation; it should be performed beforehand
%            b) returns uint32 matrix; convert it to double for further processing
%
% See also:
%    jp2_codec, TestJP2, MakePlainImage

x_size = size(Isrc,2);
y_size = size(Isrc,1);

I_cropped = int16(Isrc);

header_size_in_bytes = 0;                                                           % for large images, it is not important
target_size = round(header_size_in_bytes + (x_size * y_size * bit_depth) / (comp_ratio * 8));     % 8 for bytes instead of bits


[comp_data] = jp2_class(MakePlainImage(I_cropped), target_size, bit_depth +1);
[decoded] = jp2_class(comp_data, 0, 0);
I_new = MakePlainImage(decoded);
