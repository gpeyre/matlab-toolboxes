function c = hilbert_curve(n)

% hilbert_curve - create an hilbert curve at level n of recursion
%
%   c = hilbert_curve(n);
%
%   Copyright (c) 2004 Gabriel Peyré

if n<=0
  c = [0; 0];
else
  c0 = hilbert_curve(n-1);
  x = .5*[-1/2 + c0(2,:), c0(1,:)-1/2, c0(1,:)+1/2, -c0(2,:)+1/2];
  y = .5*[-1/2 + c0(1,:), c0(2,:)+1/2, c0(2,:)+1/2, -c0(1,:)-1/2];
  c = [x;y];
end