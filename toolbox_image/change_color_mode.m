function [M,ColorP] = change_color_mode(M,dir,options)

% change_color_mode - change of color representation
%
%   M = change_color_mode(M,dir,options);
%
%   switch between RGB and another color mode.
%
%   options.color_mode can be:
%   'rgb' (no change), 'hsv', 'ycbcr', 'pca'
%
%   For PCA you should proceed like this:
%       [Mg,options.ColorP] = change_color_mode(M0,+1,options);
%       M1 = change_color_mode(Mg,-1,options);
%
%   Copyright (c) 2007 Gabriel Peyre


options.null = 0;
color_mode = getoptions(options, 'color_mode', 'hsv');

if strcmp(color_mode,'rgb')
    return;
end

n = size(M,1);
if size(M,3)~=3
    error('Works only for color images');
end

if strcmp(color_mode, 'pca')
    nsamples = min(8000,n^2);
    if not(isfield(options, 'ColorP'))
        if dir==-1
            error('You must provide options.ColorP');
        end
        if nargout==1
            error('Output argument ColorP Missing');
        end
        % compute PCA components
        Ma = reshape(M,n^2,3);
        sel = randperm(n^2); sel = sel(1:nsamples);
        [Y,X1,v,Psi] = pca( Ma(sel,:)',3);
        ColorP = cat(2, Y, Psi(:));
    else
        ColorP = options.ColorP;
    end
    if dir==1
        M = reshape(M,n^2,3)';
        M = ColorP(:,1:3)'*(M - repmat(ColorP(:,4), 1,n^2));
        M = reshape(M',n,n,3);
    else
        M = reshape(M,n^2,3)';
        M = ColorP(:,1:3)*M + repmat(ColorP(:,4), 1,n^2);
        M = reshape(M',n,n,3);
    end
    return;
end

warning off;
if dir==1
    M = feval( ['rgb2' color_mode], M);
end
if strcmp(color_mode, 'hsv')
    M = M(:,:,3:-1:1);
end
if dir==-1
    M = feval( [color_mode '2rgb'], M);
end
warning on;





function out = rgb2ycbcr(in)
%RGB2YCBCR Convert RGB values to YCBCR color space.
%   YCBCRMAP = RGB2YCBCR(MAP) converts the RGB values in MAP to
%   the YCBCR color space.  YCBCRMAP is a M-by-3 matrix that contains
%   the YCBCR luminance (Y) and chrominance (Cb and Cr) color values as
%   columns.  Each row represents the equivalent color to the
%   corresponding row in the RGB colormap.
%
%   YCBCR = RGB2YCBCR(RGB) converts the truecolor image RGB to the
%   equivalent image in the YCBCR color space.
%
%   Class Support
%   -------------
%   If the input is an RGB image, it can be of class uint8, uint16,
%   or double; the output image is of the same class as the input 
%   image.  If the input is a colormap, the input and output colormaps 
%   are both of class double.
%
%   Examples
%   --------
%   Convert RGB image to YCbCr.
%
%      RGB = imread('board.tif');
%      YCBCR = rgb2ycbcr(RGB);
%
%   Convert RGB color space to YCbCr.
%
%      map = jet(256);
%      newmap = rgb2ycbcr(map);

%   See also NTSC2RGB, RGB2NTSC, YCBCR2RGB.

%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 1.13.4.2 $  $Date: 2003/08/23 05:54:37 $

%   Reference: 
%     C.A. Poynton, "A Technical Introduction to Digital Video", John Wiley
%     & Sons, Inc., 1996, p. 175


%initialize variables
isColormap = false;
classin = class(in);

%must reshape colormap to be m x n x 3 for transformation
if (ndims(in)==2)
  %colormap
  isColormap=true;
  colors = size(in,1);
  in = reshape(in, [colors 1 3]);
end

% set up constants for transformation
T = [65.481 128.553 24.966;...
     -37.797 -74.203 112; ...
     112 -93.786 -18.214];
offset = [16;128;128];
offset16 = 257 * offset;
fac8 = 1/255;
fac16 = 257/65535;

%initialize output
out = in;

% do transformation
switch classin
 case 'uint8'
  for p=1:3
    out(:,:,p) = imlincomb(T(p,1)*fac8,in(:,:,1),T(p,2)*fac8,in(:,:,2),...
                         T(p,3)*fac8,in(:,:,3),offset(p));
  end
 case 'uint16'
  for p=1:3
    out(:,:,p) = imlincomb(T(p,1)*fac16,in(:,:,1),T(p,2)*fac16,in(:,:,2),...
                         T(p,3)*fac16,in(:,:,3),offset16(p));
  end
 case 'double'
  % These equations transform RGB in [0,1] to YCBCR in [0, 255]
  for p=1:3
    out(:,:,p) =  T(p,1) * in(:,:,1) + T(p,2) * in(:,:,2) + T(p,3) * ...
        in(:,:,3) + offset(p);
  end
  out = out / 255;
end

if isColormap
   out = reshape(out, [colors 3 1]);
end


function rgb = ycbcr2rgb(in)
%YCBCR2RGB Convert YCbCr values to RGB color space.
%   RGBMAP = YCBCR2RGB(YCBCRMAP) converts the YCbCr values in the
%   colormap YCBCRMAP to the RGB color space. If YCBCRMAP is M-by-3 and
%   contains the YCbCr luminance (Y) and chrominance (Cb and Cr) color
%   values as columns, then RGBMAP is an M-by-3 matrix that contains
%   the red, green, and blue values equivalent to those colors.
%
%   RGB = YCBCR2RGB(YCBCR) converts the YCbCr image to the equivalent
%   truecolor image RGB.
%
%   Class Support
%   -------------
%   If the input is a YCbCr image, it can be of class uint8, uint16,
%   or double; the output image is of the same class as the input 
%   image.  If the input is a colormap, the input and output 
%   colormaps are both of class double.
%
%   Example
%   -------
%   Convert image from RGB space to YCbCr space and back.
%
%       rgb = imread('board.tif');
%       ycbcr = rgb2ycbcr(rgb);
%       rgb2 = ycbcr2rgb(ycbcr);
%
%   See also NTSC2RGB, RGB2NTSC, RGB2YCBCR.

%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 1.15.4.2 $  $Date: 2003/08/23 05:54:51 $

%   Reference:
%     Charles A. Poynton, "A Technical Introduction to Digital Video",
%     John Wiley & Sons, Inc., 1996

%initialize variables
isColormap = false;
classin = class(in);

%must reshape colormap to be m x n x 3 for transformation
if (ndims(in)==2)
  %colormap
  isColormap=true;
  colors = size(in,1);
  in = reshape(in, [colors 1 3]);
end

%initialize output
rgb = in;

% set up constants for transformation. T alone will transform YCBCR in [0,255]
% to RGB in [0,1]. We must scale T and the offsets to get RGB in the appropriate
% range for uint8 and for uint16 arrays.
T = [65.481 128.553 24.966;...
     -37.797 -74.203 112; ...
     112 -93.786 -18.214];
Tinv = T^-1;
offset = [16;128;128];
Td = 255 * Tinv;
offsetd = Tinv * offset;
T8 = Td;
offset8 = T8 * offset;
T16 = (65535/257) * Tinv;
offset16 = 65535 * Tinv * offset;

switch classin
 case 'double'
  for p = 1:3
    rgb(:,:,p) = Td(p,1) * in(:,:,1) + Td(p,2) * in(:,:,2) + ...
        Td(p,3) * in(:,:,3) - offsetd(p);
  end
  
 case 'uint8'
  for p = 1:3
    rgb(:,:,p) = imlincomb(T8(p,1),in(:,:,1),T8(p,2),in(:,:,2), ...
                           T8(p,3),in(:,:,3),-offset8(p));
  end
  
 case 'uint16'
  for p = 1:3
    rgb(:,:,p) = imlincomb(T16(p,1),in(:,:,1),T16(p,2),in(:,:,2), ...
                           T16(p,3),in(:,:,3),-offset16(p));
  end  
end

if isColormap
   rgb = reshape(rgb, [colors 3 1]);
end

if isa(rgb,'double')
  rgb = min(max(rgb,0.0),1.0);
end

%%%
%Parse Inputs
%%%
function X = parse_inputs(varargin)

checknargin(1,1,nargin,mfilename);
X = varargin{1};

if ndims(X)==2
  checkinput(X,{'uint8','uint16','double'},'real nonempty',mfilename,'MAP',1);
  if (size(X,2) ~=3 || size(X,1) < 1)
    eid = sprintf('Images:%s:invalidSizeForColormap',mfilename);
    msg = 'MAP must be a m x 3 array.';
    error(eid,'%s',msg);
  end
elseif ndims(X)==3
  checkinput(X,{'uint8','uint16','double'},'real',mfilename,'RGB',1);
  if (size(X,3) ~=3)
    eid = sprintf('Images:%s:invalidTruecolorImage',mfilename);
    msg = 'RGB must a m x n x 3 array.';
    error(eid,'%s',msg);
  end
else
  eid = sprintf('Images:%s:invalidInputSize',mfilename);
  msg = ['RGB2GRAY only accepts two-dimensional or three-dimensional ' ...
         'inputs.'];
  error(eid,'%s',msg);
end

