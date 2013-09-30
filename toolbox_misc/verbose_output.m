function verbose_output(str)

% verbose_output - display a string.
% 
%   verbose_output(str);
%
%   The string is not displayed if the global variable 'verb' is set to 0.
%   You can customize the header of the display via the global variable
%   'verb_header'.
%
%   Copyright (c) 2004 Gabriel Peyré

global verb;
global verb_header;

if isempty(verb_header)
    verb_header = '--> ';
end

if isempty(verb) || verb==1
    disp( [verb_header, str] );
end