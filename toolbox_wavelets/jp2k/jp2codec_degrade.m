function [result_dble] = jp2codec_degrade(Isrc, comp_ratio, bit_depth, Isrc_dble)
%
% Attention: a) no rounding operation; it should be performed beforehand
%            b) returns uint32 matrix; convert it to double for further
%            processing
%
% See also:
%    jp2_codec, TestJP2, MakePlainImage


%2 voies d'entree-sortie.
% voie classique :  entree de l'image dans la chaine de compression/ sortie de l'image apres degradation
% voie d'insertion dans le codeur : entree de la transformee inseree directement dans la chaine en lieu et place de
% 				   la transformee normalement calculee a partir de l'image.la sortie se fait avant
%				   la phase normale d'inversion de la transformee.
% utilit? de garder l'image en entree. L'encodeur jp2k peut se baser sur l'image originale pour definir des strategies
% optimales de codage qu'il appliquera ensuite au codage de la transformee (decoupage en tuile ...) et a son decodage.
% Dans l'utilisation qu'il est ici faite de jpeg2k, l'image d'entree n'a aucune incidence quant au codage de la transformee
% (peut etre dans une version future ?)

x_size = size(Isrc,2);
y_size = size(Isrc,1);

I_cropped = int16(Isrc);

header_size_in_bytes = 0;                                                                   % for large images, it is not important
target_size = header_size_in_bytes + (x_size * y_size * bit_depth) / (comp_ratio * 8);      % 8 for bytes instead of bits

% MakePlainImage(I_cropped)
[comp_data] = jp2_codec_encod(I_cropped, target_size, bit_depth+1, Isrc_dble);
[decoded,result_dble] = jp2_codec_encod(comp_data, 0, 0, 0);
result_dble = result_dble;


