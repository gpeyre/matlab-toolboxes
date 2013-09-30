function v = s2v(s,a)

% s2v - structure array to vector
%
%   v = s2v(s,a);
%
%   v(i) = getfield(s(i), a);
%
%   Copyright (c) 2010 Gabriel Peyre

v = zeros(length(s),1);
for i=1:length(s)
    v(i) = getfield(s(i), a);
end
