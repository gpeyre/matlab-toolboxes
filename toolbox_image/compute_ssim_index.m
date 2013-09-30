function [mssim, ssim_map] = compute_ssim_index(img1, img2, K, window, L)

%========================================================================
%SSIM Index, Version 1.0
%Copyright(c) 2003 Zhou Wang
%All Rights Reserved.
%
%The author is with Howard Hughes Medical Institute, and Laboratory
%for Computational Vision at Center for Neural Science and Courant
%Institute of Mathematical Sciences, New York University.
%
%----------------------------------------------------------------------
%Permission to use, copy, or modify this software and its documentation
%for educational and research purposes only and without fee is hereby
%granted, provided that this copyright notice and the original authors'
%names appear on all copies and supporting documentation. This program
%shall not be used, rewritten, or adapted as the basis of a commercial
%software or hardware product without first obtaining permission of the
%authors. The authors make no representations about the suitability of
%this software for any purpose. It is provided "as is" without express
%or implied warranty.
%----------------------------------------------------------------------
%
%This is an implementation of the algorithm for calculating the
%Structural SIMilarity (SSIM) index between two images. Please refer
%to the following paper:
%
%Z. Wang, A. C. Bovik, H. R. Sheikh, and E. P. Simoncelli, "Image
%quality assessment: From error visibility to structural similarity"
%IEEE Transactios on Image Processing, vol. 13, no. 4, pp.600-612,
%Apr. 2004.
%
%Kindly report any suggestions or corrections to zhouwang@ieee.org
%
%----------------------------------------------------------------------
%
%Input : (1) img1: the first image being compared
%        (2) img2: the second image being compared
%        (3) K: constants in the SSIM index formula (see the above
%            reference). defualt value: K = [0.01 0.03]
%        (4) window: local window for statistics (see the above
%            reference). default widnow is Gaussian given by
%            window = fspecial('gaussian', 11, 1.5);
%        (5) L: dynamic range of the images. default: L = 255
%
%Output: (1) mssim: the mean SSIM index value between 2 images.
%            If one of the images being compared is regarded as 
%            perfect quality, then mssim can be considered as the
%            quality measure of the other image.
%            If img1 = img2, then mssim = 1.
%        (2) ssim_map: the SSIM index map of the test image. The map
%            has a smaller size than the input images. The actual size:
%            size(img1) - size(window) + 1.
%
%Default Usage:
%   Given 2 test images img1 and img2, whose dynamic range is 0-255
%
%   [mssim ssim_map] = ssim_index(img1, img2);
%
%Advanced Usage:
%   User defined parameters. For example
%
%   K = [0.05 0.05];
%   window = ones(8);
%   L = 100;
%   [mssim ssim_map] = ssim_index(img1, img2, K, window, L);
%
%See the results:
%
%   mssim                        %Gives the mssim value
%   imshow(max(0, ssim_map).^4)  %Shows the SSIM index map
%
%========================================================================

Lmax = 1;

if (nargin < 2 | nargin > 5)
   ssim_index = -Inf;
   ssim_map = -Inf;
   return;
end

if (size(img1) ~= size(img2))
   ssim_index = -Inf;
   ssim_map = -Inf;
   return;
end

[M N] = size(img1);

if (nargin == 2)
   if ((M < 11) | (N < 11))
	   ssim_index = -Inf;
	   ssim_map = -Inf;
      return
   end
   window = fspecial('gaussian', 11, 1.5);	%
   K(1) = 0.01;								      % default settings
   K(2) = 0.03;								      %
   L = Lmax;                                  %
end

if (nargin == 3)
   if ((M < 11) | (N < 11))
	   ssim_index = -Inf;
	   ssim_map = -Inf;
      return
   end
   window = fspecial('gaussian', 11, 1.5);
   L = Lmax;
   if (length(K) == 2)
      if (K(1) < 0 | K(2) < 0)
		   ssim_index = -Inf;
   		ssim_map = -Inf;
	   	return;
      end
   else
	   ssim_index = -Inf;
   	ssim_map = -Inf;
	   return;
   end
end

if (nargin == 4)
   [H W] = size(window);
   if ((H*W) < 4 | (H > M) | (W > N))
	   ssim_index = -Inf;
	   ssim_map = -Inf;
      return
   end
   L = Lmax;
   if (length(K) == 2)
      if (K(1) < 0 | K(2) < 0)
		   ssim_index = -Inf;
   		ssim_map = -Inf;
	   	return;
      end
   else
	   ssim_index = -Inf;
   	ssim_map = -Inf;
	   return;
   end
end

if (nargin == 5)
   [H W] = size(window);
   if ((H*W) < 4 | (H > M) | (W > N))
	   ssim_index = -Inf;
	   ssim_map = -Inf;
      return
   end
   if (length(K) == 2)
      if (K(1) < 0 | K(2) < 0)
		   ssim_index = -Inf;
   		ssim_map = -Inf;
	   	return;
      end
   else
	   ssim_index = -Inf;
   	ssim_map = -Inf;
	   return;
   end
end

C1 = (K(1)*L)^2;
C2 = (K(2)*L)^2;
window = window/sum(sum(window));
img1 = double(img1);
img2 = double(img2);

mu1   = filter2(window, img1, 'valid');
mu2   = filter2(window, img2, 'valid');
mu1_sq = mu1.*mu1;
mu2_sq = mu2.*mu2;
mu1_mu2 = mu1.*mu2;
sigma1_sq = filter2(window, img1.*img1, 'valid') - mu1_sq;
sigma2_sq = filter2(window, img2.*img2, 'valid') - mu2_sq;
sigma12 = filter2(window, img1.*img2, 'valid') - mu1_mu2;

if (C1 > 0 & C2 > 0)
   ssim_map = ((2*mu1_mu2 + C1).*(2*sigma12 + C2))./((mu1_sq + mu2_sq + C1).*(sigma1_sq + sigma2_sq + C2));
else
   numerator1 = 2*mu1_mu2 + C1;
   numerator2 = 2*sigma12 + C2;
	denominator1 = mu1_sq + mu2_sq + C1;
   denominator2 = sigma1_sq + sigma2_sq + C2;
   ssim_map = ones(size(mu1));
   index = (denominator1.*denominator2 > 0);
   ssim_map(index) = (numerator1(index).*numerator2(index))./(denominator1(index).*denominator2(index));
   index = (denominator1 ~= 0) & (denominator2 == 0);
   ssim_map(index) = numerator1(index)./denominator1(index);
end

mssim = mean(ssim_map(:));

return



function h = fspecial(varargin)
%FSPECIAL Create 2-D special filters.
%   H = FSPECIAL(TYPE) creates a two-dimensional filter H of the
%   specified type. Possible values for TYPE are:
%
%     'average'   averaging filter
%     'disk'      circular averaging filter
%     'gaussian'  Gaussian lowpass filter
%     'laplacian' filter approximating the 2-D Laplacian operator
%     'log'       Laplacian of Gaussian filter
%     'motion'    motion filter
%     'prewitt'   Prewitt horizontal edge-emphasizing filter
%     'sobel'     Sobel horizontal edge-emphasizing filter
%     'unsharp'   unsharp contrast enhancement filter
%
%   Depending on TYPE, FSPECIAL may take additional parameters
%   which you can supply.  These parameters all have default
%   values. 
%
%   H = FSPECIAL('average',HSIZE) returns an averaging filter H of size
%   HSIZE. HSIZE can be a vector specifying the number of rows and columns in
%   H or a scalar, in which case H is a square matrix.
%   The default HSIZE is [3 3].
%
%   H = FSPECIAL('disk',RADIUS) returns a circular averaging filter
%   (pillbox) within the square matrix of side 2*RADIUS+1.
%   The default RADIUS is 5.
%
%   H = FSPECIAL('gaussian',HSIZE,SIGMA) returns a rotationally
%   symmetric Gaussian lowpass filter  of size HSIZE with standard
%   deviation SIGMA (positive). HSIZE can be a vector specifying the
%   number of rows and columns in H or a scalar, in which case H is a
%   square matrix.
%   The default HSIZE is [3 3], the default SIGMA is 0.5.
%
%   H = FSPECIAL('laplacian',ALPHA) returns a 3-by-3 filter
%   approximating the shape of the two-dimensional Laplacian
%   operator. The parameter ALPHA controls the shape of the
%   Laplacian and must be in the range 0.0 to 1.0.
%   The default ALPHA is 0.2.
%
%   H = FSPECIAL('log',HSIZE,SIGMA) returns a rotationally symmetric
%   Laplacian of Gaussian filter of size HSIZE with standard deviation
%   SIGMA (positive). HSIZE can be a vector specifying the number of rows
%   and columns in H or a scalar, in which case H is a square matrix.
%   The default HSIZE is [5 5], the default SIGMA is 0.5.
%
%   H = FSPECIAL('motion',LEN,THETA) returns a filter to approximate, once
%   convolved with an image, the linear motion of a camera by LEN pixels,
%   with an angle of THETA degrees in a counter-clockwise direction. The
%   filter becomes a vector for horizontal and vertical motions.  The
%   default LEN is 9, the default THETA is 0, which corresponds to a
%   horizontal motion of 9 pixels.
%
%   H = FSPECIAL('prewitt') returns 3-by-3 filter that emphasizes
%   horizontal edges by approximating a vertical gradient. If you need to
%   emphasize vertical edges, transpose the filter H: H'.
%
%       [1 1 1;0 0 0;-1 -1 -1].
%
%   H = FSPECIAL('sobel') returns 3-by-3 filter that emphasizes
%   horizontal edges utilizing the smoothing effect by approximating a
%   vertical gradient. If you need to emphasize vertical edges, transpose
%   the filter H: H'.
%
%       [1 2 1;0 0 0;-1 -2 -1].
%
%   H = FSPECIAL('unsharp',ALPHA) returns a 3-by-3 unsharp contrast
%   enhancement filter. FSPECIAL creates the unsharp filter from the
%   negative of the Laplacian filter with parameter ALPHA. ALPHA controls
%   the shape of the Laplacian and must be in the range 0.0 to 1.0.
%   The default ALPHA is 0.2.
%
%   Class Support
%   -------------
%   H is of class double.
%
%   Example
%   -------
%      I = imread('cameraman.tif');
%      subplot(2,2,1);imshow(I);title('Original Image'); 
%      H = fspecial('motion',20,45);
%      MotionBlur = imfilter(I,H,'replicate');
%      subplot(2,2,2);imshow(MotionBlur);title('Motion Blurred Image');
%      H = fspecial('disk',10);
%      blurred = imfilter(I,H,'replicate');
%      subplot(2,2,3);imshow(blurred);title('Blurred Image');
%      H = fspecial('unsharp');
%      sharpened = imfilter(I,H,'replicate');
%      subplot(2,2,4);imshow(sharpened);title('Sharpened Image');
%       
%   See also CONV2, EDGE, FILTER2, FSAMP2, FWIND1, FWIND2, IMFILTER.

%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 5.28.4.2 $  $Date: 2003/01/26 05:55:24 $

[type, p2, p3] = ParseInputs(varargin{:});

switch type
  case 'average' % Smoothing filter
     siz = p2;
     h   = ones(siz)/prod(siz);

  case 'disk' % Disk filter
     rad   = p2;
     crad  = ceil(rad-0.5);
     [x,y] = meshgrid(-crad:crad,-crad:crad);
     maxxy = max(abs(x),abs(y));
     minxy = min(abs(x),abs(y));
     m1 = (rad^2 <  (maxxy+0.5).^2 + (minxy-0.5).^2).*(minxy-0.5) + ...
          (rad^2 >= (maxxy+0.5).^2 + (minxy-0.5).^2).* ...
	        sqrt(rad^2 - (maxxy + 0.5).^2);
     m2 = (rad^2 >  (maxxy-0.5).^2 + (minxy+0.5).^2).*(minxy+0.5) + ...
          (rad^2 <= (maxxy-0.5).^2 + (minxy+0.5).^2).* ...
           sqrt(rad^2 - (maxxy - 0.5).^2);
     sgrid = (rad^2*(0.5*(asin(m2/rad) - asin(m1/rad)) + ...
             0.25*(sin(2*asin(m2/rad)) - sin(2*asin(m1/rad)))) - ...
             (maxxy-0.5).*(m2-m1) + (m1-minxy+0.5)) ... 
	          .*((((rad^2 < (maxxy+0.5).^2 + (minxy+0.5).^2) & ...
             (rad^2 > (maxxy-0.5).^2 + (minxy-0.5).^2)) | ...
	          ((minxy==0)&(maxxy-0.5 < rad)&(maxxy+0.5>=rad))));
     sgrid = sgrid + ((maxxy+0.5).^2 + (minxy+0.5).^2 < rad^2);
     sgrid(crad+1,crad+1) = min(pi*rad^2,pi/2);
     if ((crad>0) & (rad > crad-0.5) & (rad^2 < (crad-0.5)^2+0.25)) 
        m1  = sqrt(rad^2 - (crad - 0.5).^2);
	     m1n = m1/rad;
        sg0 = 2*(rad^2*(0.5*asin(m1n) + 0.25*sin(2*asin(m1n)))-m1*(crad-0.5));
        sgrid(2*crad+1,crad+1) = sg0;
        sgrid(crad+1,2*crad+1) = sg0;
        sgrid(crad+1,1)        = sg0;
        sgrid(1,crad+1)        = sg0;
        sgrid(2*crad,crad+1)   = sgrid(2*crad,crad+1) - sg0;
        sgrid(crad+1,2*crad)   = sgrid(crad+1,2*crad) - sg0;
        sgrid(crad+1,2)        = sgrid(crad+1,2)      - sg0;
        sgrid(2,crad+1)        = sgrid(2,crad+1)      - sg0;
     end
     sgrid(crad+1,crad+1) = min(sgrid(crad+1,crad+1),1);
     h = sgrid/sum(sgrid(:));

  case 'gaussian' % Gaussian filter

     siz   = (p2-1)/2;
     std   = p3;
     
     [x,y] = meshgrid(-siz(2):siz(2),-siz(1):siz(1));
     arg   = -(x.*x + y.*y)/(2*std*std);

     h     = exp(arg);
     h(h<eps*max(h(:))) = 0;

     sumh = sum(h(:));
     if sumh ~= 0,
       h  = h/sumh;
     end;
     
  case 'laplacian' % Laplacian filter
     alpha = p2;
     alpha = max(0,min(alpha,1));
     h1    = alpha/(alpha+1); h2 = (1-alpha)/(alpha+1);
     h     = [h1 h2 h1;h2 -4/(alpha+1) h2;h1 h2 h1];

  case 'log' % Laplacian of Gaussian
     % first calculate Gaussian
     siz   = (p2-1)/2;
     std2   = p3^2;
     
     [x,y] = meshgrid(-siz(2):siz(2),-siz(1):siz(1));
     arg   = -(x.*x + y.*y)/(2*std2);
     
     h     = exp(arg);
     h(h<eps*max(h(:))) = 0;

     sumh = sum(h(:));
     if sumh ~= 0,
       h  = h/sumh;
     end;
     % now calculate Laplacian     
     h1 = h.*(x.*x + y.*y - 2*std2)/(std2^2);
     h     = h1 - sum(h1(:))/prod(p2); % make the filter sum to zero
  
  case 'motion' % Motion filter uses bilinear interpolation
     len = max(1,p2);
     half = (len-1)/2;% rotate half length around center
     phi = mod(p3,180)/180*pi;

     cosphi = cos(phi);
     sinphi = sin(phi);
     xsign = sign(cosphi);
     linewdt = 1;

     % define mesh for the half matrix, eps takes care of the right size
     % for 0 & 90 rotation
     sx = fix(half*cosphi + linewdt*xsign - len*eps);
     sy = fix(half*sinphi + linewdt - len*eps);
     [x y] = meshgrid([0:xsign:sx],[0:sy]);

     % define shortest distance from a pixel to the rotated line 
     dist2line = (y*cosphi-x*sinphi);% distance perpendicular to the line

     rad = sqrt(x.^2 + y.^2);
     % find points beyond the line's end-point but within the line width
     lastpix = find((rad >= half)&(abs(dist2line)<=linewdt));
     %distance to the line's end-point parallel to the line 
     x2lastpix = half - abs((x(lastpix) + dist2line(lastpix)*sinphi)/cosphi);

     dist2line(lastpix) = sqrt(dist2line(lastpix).^2 + x2lastpix.^2);
     dist2line = linewdt + eps - abs(dist2line);
     dist2line(dist2line<0) = 0;% zero out anything beyond line width

     % unfold half-matrix to the full size
     h = rot90(dist2line,2);
     h(end+[1:end]-1,end+[1:end]-1) = dist2line;
     h = h./(sum(h(:)) + eps*len*len);

     if cosphi>0,
       h = flipud(h);
     end
     
  case 'prewitt' % Prewitt filter
     h = [1 1 1;0 0 0;-1 -1 -1];

  case 'sobel' % Sobel filter
     h = [1 2 1;0 0 0;-1 -2 -1];

  case 'unsharp' % Unsharp filter
     alpha = p2;
     h     = [0 0 0;0 1 0;0 0 0] - fspecial('laplacian',alpha);

  end
  

%%%
%%% ParseInputs
%%%
function [type, p2, p3] = ParseInputs(varargin)

% default values
type      = '';
p2        = [];
p3        = [];

% Check the number of input arguments.
% checknargin(1,3,nargin,mfilename);

% Determine filter type from the user supplied string.
type = varargin{1};
% type = checkstrs(type,{'gaussian','sobel','prewitt','laplacian','log',...
%                       'average','unsharp','disk','motion'},mfilename,'TYPE',1);
  
% default values
switch type
	case 'average'
      p2 = [3 3];  % siz
      
   case 'disk'
      p2 = 5;      % rad
      
   case 'gaussian'
      p2 = [3 3];  % siz
      p3 = 0.5;    % std
      
   case {'laplacian', 'unsharp'}
      p2 = 1/5;    % alpha
      
   case 'log'
      p2 = [5 5];  % siz
      p3 = 0.5;    % std
      
   case 'motion'
      p2 = 9;     % len
      p3 = 0;      % theta
   end
   

switch nargin
    case 1
        % FSPECIAL('average')
        % FSPECIAL('disk')
        % FSPECIAL('gaussian')
        % FSPECIAL('laplacian')
        % FSPECIAL('log')
        % FSPECIAL('motion')
        % FSPECIAL('prewitt')
        % FSPECIAL('sobel')
        % FSPECIAL('unsharp')
        % Nothing to do here; the default values have 
        % already been assigned.        
        
    case 2
       % FSPECIAL('average',N)
       % FSPECIAL('disk',RADIUS)
       % FSPECIAL('gaussian',N)
       % FSPECIAL('laplacian',ALPHA)
       % FSPECIAL('log',N)
       % FSPECIAL('motion',LEN)
       % FSPECIAL('unsharp',ALPHA)
       p2 = varargin{2};
 
       switch type
          case {'sobel','prewitt'}
              msg = sprintf('%s: Too many arguments for this type of filter.', upper(mfilename));
              eid = sprintf('Images:%s:tooManyArgsForThisFilter', mfilename);
              error(eid,msg);
          case {'laplacian','unsharp'}
              checkinput(p2,{'double'},{'nonnegative','real',...
                         'nonempty','finite','scalar'},...
                         mfilename,'ALPHA',2);
              if  p2 > 1
                  msg = sprintf('%s: ALPHA should be less than or equal 1 and greater than 0.', upper(mfilename));
                  eid = sprintf('Images:%s:outOfRangeAlpha', mfilename);
                  error(eid,msg);
              end
          case {'disk','motion'}
              checkinput(p2,{'double'},{'positive','finite','real','nonempty','scalar'},mfilename,'RADIUS or LEN',2);
          case {'gaussian','log','average'}
              checkinput(p2,{'double'},{'positive','finite','real','nonempty','integer'},mfilename,'HSIZE',2);
              if prod(size(p2)) > 2
                  msg = 'HSIZE should have 1 or 2 elements.';
                  eid = sprintf('Images:%s:wrongSizeN', mfilename);
                  error(eid,msg);
              elseif (prod(size(p2))==1)
                  p2 = [p2 p2]; 
              end
       end       

       
    case 3
       % FSPECIAL('gaussian',N,SIGMA)
       % FSPECIAL('log',N,SIGMA)
       % FSPECIAL('motion',LEN,THETA)
       p2 = varargin{2};
       p3 = varargin{3};
       
       switch type
          case 'motion'
%              checkinput(p2,{'double'},{'positive','finite','real','nonempty','scalar'},mfilename,'LEN',2);
%              checkinput(p3,{'double'},{'real','nonempty','finite','scalar'},mfilename,'THETA',3);
          case {'gaussian','log'}
%              checkinput(p2,{'double'},{'positive','finite','real','nonempty','integer'},mfilename,'N',2);
%              checkinput(p3,{'double'},{'positive','finite','real','nonempty','scalar'},mfilename,'SIGMA',3);
              if prod(size(p2)) > 2
                  msg = sprintf('%s: size(N) should be less than or equal 2.', upper(mfilename));
                  eid = sprintf('Images:%s:wrongSizeN', mfilename);
                  error(eid,msg);
              elseif (prod(size(p2))==1)
                  p2 = [p2 p2]; 
              end
          otherwise   
              msg = sprintf('%s: Too many arguments for this type of filter.', upper(mfilename));
              eid = sprintf('Images:%s:tooManyArgsForThisFilter', mfilename);
              error(eid,msg);
      end
end
