function b = imnoise(varargin)
%IMNOISE Add noise to image.
%   J = IMNOISE(I,TYPE,...) Add noise of a given TYPE to the intensity image
%   I. TYPE is a string that can have one of these values:
%
%       'gaussian'       Gaussian white noise with constant
%                        mean and variance
%
%       'localvar'       Zero-mean Gaussian white noise 
%                        with an intensity-dependent variance
%
%       'poisson'        Poisson noise
%
%       'salt & pepper'  "On and Off" pixels
%
%       'speckle'        Multiplicative noise
%
%   Depending on TYPE, you can specify additional parameters to IMNOISE. All
%   numerical parameters are normalized; they correspond to operations with
%   images with intensities ranging from 0 to 1.
%   
%   J = IMNOISE(I,'gaussian',M,V) adds Gaussian white noise of mean M and
%   variance V to the image I. When unspecified, M and V default to 0 and
%   0.01 respectively.
%   
%   J = imnoise(I,'localvar',V) adds zero-mean, Gaussian white noise of
%   local variance, V, to the image I.  V is an array of the same size as I.
%
%   J = imnoise(I,'localvar',IMAGE_INTENSITY,VAR) adds zero-mean, Gaussian
%   noise to an image, I, where the local variance of the noise is a
%   function of the image intensity values in I.  IMAGE_INTENSITY and VAR
%   are vectors of the same size, and PLOT(IMAGE_INTENSITY,VAR) plots the
%   functional relationship between noise variance and image intensity.
%   IMAGE_INTENSITY must contain normalized intensity values ranging from 0
%   to 1.
%
%   J = IMNOISE(I,'poisson') generates Poisson noise from the data instead
%   of adding artificial noise to the data. In order to respect Poisson
%   statistics, the intensities of uint8 and uint16 images must correspond
%   to the number of photons (or any other quanta of information).
%   Double-precision images are used when the number of photons per pixel
%   can be much larger than 65535 (but less than 10^12); the intensities
%   values vary between 0 and 1 and correspond to the number of photons
%   divided by 10^12.
%
%   J = IMNOISE(I,'salt & pepper',D) adds "salt and pepper" noise to the
%   image I, where D is the noise density.  This affects approximately
%   D*numel(I) pixels. The default for D is 0.05.
%   
%   J = IMNOISE(I,'speckle',V) adds multiplicative noise to the image I,
%   using the equation J = I + n*I, where n is uniformly distributed random
%   noise with mean 0 and variance V. The default for V is 0.04.
%
%   Note
%   ----
%   The mean and variance parameters for 'gaussian', 'localvar', and
%   'speckle' noise types are always specified as if for a double image
%   in the range [0, 1].  If the input image is of class uint8 or uint16,
%   the imnoise function converts the image to double, adds noise
%   according to the specified type and parameters, and then converts the
%   noisy image back to the same class as the input.
%   
%   Class Support
%   -------------
%   I can be of class uint8, uint16, or double.Ê The output image J is of
%   the same class as I.  If I has more than two dimensions it is treated
%   as a multidimensional intensity image and not as an RGB image.
%
%   Example
%   -------
%        I = imread('eight.tif');
%        J = imnoise(I,'salt & pepper', 0.02);
%        imview(I), imview(J)
%
%   See also IMNSTATS, RAND, RANDN.

%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 5.20.4.3 $  $Date: 2003/08/01 18:09:05 $

  [a, code, classIn, classChanged, p3, p4] = ParseInputs(varargin{:});

  clear varargin;
  sizeA = size(a);

  switch code
   case 'gaussian' % Gaussian white noise
    b = a + sqrt(p4)*randn(sizeA) + p3;
    
   case 'localvar_1' % Gaussian white noise with variance varying locally
                     % imnoise(a,'localvar',v)
                     % v is local variance array
    b = a + sqrt(p3).*randn(sizeA); % Use a local variance array
    
   case 'localvar_2' % Gaussian white noise with variance varying locally
                     % Use an empirical intensity-variance relation
    intensity = p3(:); % Use an empirical intensity-variance relation
    var       = p4(:);
    minI  = min(intensity);
    maxI  = max(intensity);
    b     = min(max(a,minI),maxI);
    b     = reshape(interp1(intensity,var,b(:)),sizeA);
    b     = a + sqrt(b).*randn(sizeA);
    
   case 'poisson' % Poisson noise
    switch classIn
     case 'uint8'
      a = round(a*255); 
     case 'uint16'
      a = round(a*65535);
     case 'double'
      a = round(a*10^12); % Recalibration
    end
    
    a = a(:);

    %  (Monte-Carlo Rejection Method) Ref. Numerical 
    %  Recipes in C, 2nd Edition, Press, Teukolsky, 
    %  Vetterling, Flannery (Cambridge Press)
    
    b=zeros(size(a));
    idx1=find(a<50); % Cases where pixel intensities are less than 50 units
    if (~isempty(idx1))
      g=exp(-a(idx1));
      em=-ones(size(g));
      t=ones(size(g));
      idx2= (1:length(idx1))';
      while ~isempty(idx2)
        em(idx2)=em(idx2)+1;
        t(idx2)=t(idx2).*rand(size(idx2));
        idx2 = idx2(t(idx2) > g(idx2));
      end
      b(idx1)=em;
    end

    % For large pixel intensities the Poisson pdf becomes 
    % very similar to a Gaussian pdf of mean and of variance
    % equal to the local pixel intensities. Ref. Mathematical Methods
    % of Physics, 2nd Edition, Mathews, Walker (Addison Wesley)
    idx1=find(a>=50); % Cases where pixel intensities are more than 49 units
    if (~isempty(idx1))
      b(idx1)=round(a(idx1)+sqrt(a(idx1)).*randn(size(idx1)));
    end
    
    b = reshape(b,sizeA);
    
   case 'salt & pepper' % Salt & pepper noise
    b = a;
    x = rand(sizeA);
    d = find(x < p3/2);
    b(d) = 0; % Minimum value
    d = find(x >= p3/2 & x < p3);
    b(d) = 1; % Maximum (saturated) value
    
   case 'speckle' % Speckle (multiplicative) noise
    b = a + sqrt(12*p3)*a.*(rand(sizeA)-.5);
    
  end

  % Truncate the output array data if necessary
  if strcmp(code,{'poisson'})
    switch classIn
     case 'uint8'
      b = uint8(b); 
     case 'uint16'
      b = uint16(b);
     case 'double'
      b = min(b/10^12,1); 
    end
  else    
    b = max(0,min(b,1));
    % The output class should be the same as the input class
    if classChanged,
      b = changeclass(classIn, b);
    end
  end


%%%
%%% ParseInputs
%%%
function [a, code, classIn, classChanged, p3, p4, msg] = ParseInputs(varargin)

% Initialization
a            = [];
code         = 'gaussian';
classIn      = '';
classChanged = [];
p3           = [];
p4           = [];
msg = '';

% Check the number of input arguments.

% checknargin(1,4,nargin,mfilename);

% Check the input-array type.
a = varargin{1};
% checkinput(a, {'uint8','uint16','double'}, '', mfilename, 'I', 1);

% Change class to double
classIn = class(a);
classChanged = 0;
if ~isa(a, 'double')
  a = im2double(a);
  classChanged = 1;
else
  % Check for valid image I
  if ~(isNonnegativeReal(a(:)) && all(a(:)<=1))
    eid = sprintf('Images:%s:invalidDoubleImage',mfilename);
    msg = 'A valid double image must have values between zero and one.';
    error(eid,'%s',msg);
  end
end

% Check the noise type.
if nargin > 1
  if ~ischar(varargin{2})
    eid = sprintf('Images:%s:invalidNoiseType',mfilename);
    msg = 'TYPE must be a character string.';
    error(eid,'%s',msg);
  end
  
  % Preprocess noise type string to detect abbreviations.
  allStrings = {'gaussian', 'salt & pepper', 'speckle',...
                'poisson','localvar'};
  idx = strmatch(lower(varargin{2}), allStrings);
  switch length(idx)
   case 0
    eid = sprintf('Images:%s:unknownNoiseType',mfilename);
    msg = sprintf('Unknown noise type: ''%s''.', varargin{2});
    error(eid,'%s',msg);
   case 1
    code = allStrings{idx};
   otherwise
    eid = sprintf('Images:%s:ambiguousNoiseType',mfilename); 
    msg = sprintf('Ambiguous noise type: ''%s''.', varargin{2});
    error(eid,'%s',msg); 
  end
else
  code = 'gaussian';  % default noise type
end 

switch code
 case 'poisson'
  if nargin > 2
    eid = sprintf('Images:%s:tooManyPoissonInputs',mfilename);
    msg = 'Too many inputs for ''poisson'' noise.';
    error(eid,'%s',msg);
  end
  
 case 'gaussian'
  p3 = 0;     % default mean
  p4 = 0.01;  % default variance
  
  if nargin > 2
    p3 = varargin{3};
    if ~isRealScalar(p3)
      eid = sprintf('Images:%s:invalidMean',mfilename);
      msg = 'For ''gaussian'' noise, M must be a real scalar.';
      error(eid,'%s',msg);
    end
  end
  
  if nargin > 3
    p4 = varargin{4};
    if ~isNonnegativeRealScalar(p4)
      eid = sprintf('Images:%s:invalidVariance',mfilename);      
      msg = 'For ''gaussian'' noise, V must be a real nonnegative scalar.';
      error(eid,'%s',msg); 
    end
  end
  
  if nargin > 4
    eid = sprintf('Images:%s:tooManyGaussianInputs',mfilename);      
    msg = 'Too many inputs for ''gaussian'' noise.';
    error(eid,'%s',msg); 
    return;
  end
  
 case 'salt & pepper'
  p3 = 0.05;   % default density
  
  if nargin > 2
    p3 = varargin{3};
    if ~isNonnegativeRealScalar(p3) || (p3 > 1)
      eid = sprintf('Images:%s:invalidNoiseDensity',mfilename);            
      msg1 = 'For ''salt & pepper'' noise, D must be a real nonnegative ';
      msg2 = 'scalar less than or equal to 1.';
      msg = sprintf('%s\n%s',msg1,msg2);
      error(eid,'%s',msg);
    end
    
    if nargin > 3
      eid = sprintf('Images:%s:tooManySaltAndPepperInputs',mfilename);      
      msg = 'Too many inputs for ''salt & pepper'' noise.';
      error(eid,'%s',msg);
    end
  end
  
 case 'speckle'
  p3 = 0.05;    % default variance
  
  if nargin > 2
    p3 = varargin{3};
    if ~isNonnegativeRealScalar(p3)
      eid = sprintf('Images:%s:invalidVariance',mfilename);
      msg = 'For ''speckle'' noise, V must be a nonnegative real scalar.';
      error(eid,'%s',msg);
    end
  end
  
  if nargin > 3
    eid = sprintf('Images:%s:tooManySpeckleInputs',mfilename);
    msg = 'Too many inputs for ''speckle'' noise.';
    error(eid,'%s',msg);
  end
  
 case 'localvar'
  if nargin < 3
    eid = sprintf('Images:%s:toofewLocalVarInputs',mfilename);
    msg = 'Too few inputs for ''localvar'' noise.';
    error(eid,'%s',msg);
    
  elseif nargin == 3
    % IMNOISE(a,'localvar',v)
    code = 'localvar_1';
    p3 = varargin{3};
    if ~isNonnegativeReal(p3) || ~isequal(size(p3),size(a))
      eid = sprintf('Images:%s:invalidLocalVariance',mfilename);
      msg1 = 'For the ''localvar'' noise syntax, V must contain ';
      msg2 = 'only nonnegative real values and be the same size as A.';
      msg = sprintf('%s\n%s',msg1,msg2);
      error(eid,'%s',msg);
    end
    
  elseif nargin == 4
    % IMNOISE(a,'localvar',IMAGE_INTENSITY,NOISE_VARIANCE)
    code = 'localvar_2';
    p3 = varargin{3};
    p4 = varargin{4};
    
    if ~isNonnegativeRealVector(p3) || (any(p3) > 1)
      eid = sprintf('Images:%s:invalidImageIntensity',mfilename);
      msg1 = 'For ''localvar'' noise, IMAGE_INTENSITY must be a ';
      msg2 = 'nonnegative real vector less than or equal to 1.';
      msg = sprintf('%s\n%s',msg1,msg2);
      error(eid,'%s',msg);
    end
    
    if ~isNonnegativeRealVector(p4)
      eid = sprintf('Images:%s:invalidLocalVariance',mfilename);
      msg1 = 'For ''localvar'' noise, NOISE_VARIANCE must be a ';
      msg2 = 'nonnegative real vector.';
      msg = sprintf('%s\n%s',msg1,msg2);
      error(eid,'%s',msg);
    end
    
    if ~isequal(size(p3),size(p4))
      eid = sprintf('Images:%s:invalidSize',mfilename);
      msg1 = 'For ''localvar'' noise, IMAGE_INTENSITY and '; 
      msg2 = 'NOISE_VARIANCE must be the same size.';
      msg = sprintf('%s\n%s',msg1,msg2);
      error(eid,'%s',msg);
    end
    
  else
    eid = sprintf('Images:%s:tooManyLocalVarInputs',mfilename);
    msg = 'Too many inputs for ''localvar'' noise.';
    error(eid,'%s',msg);
  end
  
end

%%%
%%% isReal
%%%
function t = isReal(P)
%   isReal(P) returns 1 if P contains only real  
%   numbers and returns 0 otherwise.
%
  isFinite  = all(isfinite(P(:)));
  t = isreal(P) && isFinite && ~isempty(P);


%%%
%%% isNonnegativeReal
%%%
function t = isNonnegativeReal(P)
%   isNonnegativeReal(P) returns 1 if P contains only real  
%   numbers greater than or equal to 0 and returns 0 otherwise.
%
  t = isReal(P) && all(P(:)>=0);


%%%
%%% isRealScalar
%%%
function t = isRealScalar(P)
%   isRealScalar(P) returns 1 if P is a real, 
%   scalar number and returns 0 otherwise.
%
  t = isReal(P) && (numel(P)==1);


%%%
%%% isNonnegativeRealScalar
%%%
function t = isNonnegativeRealScalar(P)
%   isNonnegativeRealScalar(P) returns 1 if P is a real, 
%   scalar number greater than 0 and returns 0 otherwise.
%
  t = isReal(P) && all(P(:)>=0) && (numel(P)==1);


%%%
%%% isVector
%%%
function t = isVector(P)
%   isVector(P) returns 1 if P is a vector and returns 0 otherwise.
%
  t = ((numel(P) >= 2) && ((size(P,1) == 1) || (size(P,2) == 1)));


%%%
%%% isNonnegativeRealVector
%%%
function t = isNonnegativeRealVector(P)
%   isNonnegativeRealVector(P) returns 1 if P is a real, 
%   vector greater than 0 and returns 0 otherwise.
%
  t = isReal(P) && all(P(:)>=0) && isVector(P);
