function [y,I,J] = remove_doublon(x)

% remove_doublon - remove doublon from a vector (synonym with UNIQUE).
%
%   [y,I,J] = remove_doublon(x);
%
%   Copyright (c) 2004 Gabriel Peyré

[y,I,J] = unique(x);
return;

%% OLD CODE
x = reshape(x,1,prod(size(x)));
y = [];
I = [];
k = 0;
for i=x
    k = k+1;
    if isempty(find(y==i))
        y = [y,i];
        I = [I,k];
    end
end