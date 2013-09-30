function [out_data, out_data_ex] = jp2_codec(in_data, compressed_size, bits_per_pixel, coeff_amp)

% this function calls via the compiled function jp2_codec.mexglx
% one of the subroutines C from the file liw_jp2_dll.c:
%   liw_jp2_encode_ex
%   liw_jp2_decode_ex
%   liw_jp2_info_ex
%
% In order to compile this function, one should type
% > compile_jp2_codec
%
% The generated Matlab function ptortotype is:
%
% function [out_data, out_data_ex] = jp2_codec(in_data, compressed_size, bits_per_pixel, coeff_amp)
%
% Inputs:
% ------
%   in_data : a color image to be compressed        ( 3xXxY uint8 matrix )
%          or a grayscale image to be compressed    ( XxY int8, int16, int32  matrix )
%          or a compressed code to be decoded       ( 1xN or Nx1 uint8 vector )
%
%   Attention: grayscale image must be transposed before and after calling this function
%              Use MakePlainImage function to do the same operation for colour images
%
%   compressed_size: 0 if decoding
%                    desired size of compressed data, in bytes, if encoding
%
%   bits_per_pixel:  ignored if decoding
%                    maximum bits per pixel per component
%
%   coeff_amp:     coefficient amplification table
%
%   out_data_ex:    some debugging JPEG-2000 information, do not use.
%
% Outputs:
% -------
%   out_data : an image or a compressed code (in the same format as in_data)
%
% See also:
% ---------
%   MakePlainImage
%   jp2_degrade
%
