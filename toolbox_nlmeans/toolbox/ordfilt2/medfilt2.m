function b = medfilt2(varargin)
%MEDFILT2 Perform 2-D median filtering.
%   B = MEDFILT2(A,[M N]) performs median filtering of the matrix
%   A in two dimensions. Each output pixel contains the median
%   value in the M-by-N neighborhood around the corresponding
%   pixel in the input image. MEDFILT2 pads the image with zeros
%   on the edges, so the median values for the points within 
%   [M N]/2 of the edges may appear distorted.
%
%   B = MEDFILT2(A) performs median filtering of the matrix A
%   using the default 3-by-3 neighborhood.
%
%   B = MEDFILT2(...,PADOPT) controls how the matrix boundaries
%   are padded.  PADOPT may be 'zeros' (the default),
%   'symmetric', or 'indexed'. If PADOPT is 'zeros', A is padded
%   with zeros at the boundaries. If PADOPT is 'symmetric', A is
%   symmetrically extended at the boundaries. If PADOPT is
%   'indexed', A is padded with ones if it is double; otherwise
%   it is padded with zeros.
%
%   Class Support
%   -------------
%   The input image A can be logical or numeric (unless the 
%   'indexed' syntax is used, in which case A cannot be of class 
%   uint16).  The output image B is of the same class as A.
%
%   Remarks
%   -------
%   If the input image A is of integer class, all of the output
%   values are returned as integers. If the number of
%   pixels in the neighborhood (i.e., M*N) is even, some of the
%   median values may not be integers. In these cases, the
%   fractional parts are discarded. Logical input is treated
%   similarly.
%
%   Example
%   -------
%       I = imread('eight.tif');
%       J = imnoise(I,'salt & pepper',0.02);
%       K = medfilt2(J);
%       imview(J), imview(K)
%
%   See also FILTER2, ORDFILT2, WIENER2.

%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 5.18.4.6 $  $Date: 2003/08/23 05:53:02 $

[a, mn, padopt] = parse_inputs(varargin{:});

domain = ones(mn);
if (rem(prod(mn), 2) == 1)
    order = (prod(mn)+1)/2;
    b = ordfilt2(a, order, domain, padopt);
else
    order1 = prod(mn)/2;
    order2 = order1+1;
    b = ordfilt2(a, order1, domain, padopt);
    b2 = ordfilt2(a, order2, domain, padopt);
	if islogical(b)
		b = b | b2;
	else
		b =	imlincomb(0.5, b, 0.5, b2);
	end
end


%%%
%%% Function parse_inputs
%%%
function [a, mn, padopt] = parse_inputs(varargin)
% checknargin(1,4,nargin,mfilename);

% There are several grandfathered syntaxes we have to catch
% and parse successfully, so we're going to use a strategy
% that's a little different that usual.
%
% First, scan the input argument list for strings.  The
% string 'indexed', 'zeros', or 'symmetric' can appear basically
% anywhere after the first argument.
%
% Second, delete the strings from the argument list.
%
% The remaining argument list can be one of the following:
% MEDFILT2(A)
% MEDFILT2(A,[M N])
% MEDFILT2(A,[M N],[Mb Nb])
%
% Any syntax in which 'indexed' is followed by other arguments
% is grandfathered.  Any syntax in which [Mb Nb] appears is
% grandfathered.
%
% -sle, March 1998

a = varargin{1};

charLocation = [];
for k = 2:nargin
    if (ischar(varargin{k}))
        charLocation = [charLocation k];
    end
end

if (length(charLocation) > 1)
    % More than one string in input list
    eid = 'Images:medfilt2:tooManyStringInputs';
    error(eid,'%s','Too many input string arguments.');
elseif isempty(charLocation)
    % No string specified
    padopt = 'zeros';
else
    options = {'indexed', 'zeros', 'symmetric'};

    padopt = checkstrs(varargin{charLocation}, options, mfilename, ...
                       'PADOPT', charLocation);
    
    varargin(charLocation) = [];
end

if (strcmp(padopt, 'indexed'))
    if (isa(a,'double'))
        padopt = 'ones';
    else
        padopt = 'zeros';
    end
end

if length(varargin) == 1,
  mn = [3 3];% default
elseif length(varargin) >= 2,
  mn = varargin{2}(:).';
  if size(mn,2)~=2,
    msg = 'MEDFILT2(A,[M N]): Second argument must consist of two integers.';
    eid = 'Images:medfilt2:secondArgMustConsistOfTwoInts';
    error(eid, msg);
  elseif length(varargin) > 2,
    msg = ['MEDFILT2(A,[M N],[Mb Nb],...) is an obsolete syntax. [Mb Nb]' ...
             ' argument is ignored.'];
    wid = 'Images:medfilt2:obsoleteSyntax';
    warning(wid, msg);
  end
end

% The grandfathered [Mb Nb] argument, if present, is ignored.
