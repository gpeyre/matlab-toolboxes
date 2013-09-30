function out = mkconstarray(class, value, size)
%MKCONSTARRAY creates a constant array of a specified numeric class.
%   A = MKCONSTARRAY(CLASS, VALUE, SIZE) creates a constant array 
%   of value VALUE and of size SIZE.

%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 1.8.4.1 $  $Date: 2003/01/26 06:00:35 $

out = repmat(feval(class, value), size);

