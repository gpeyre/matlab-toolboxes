function imagesc_log(M, lambda, x, y)

[n,p] = size(M);

if nargin<2
    lambda = 2;
end
if nargin<3
    x = 1:n;
end
if nargin<4
    y = 1:p;
end

m = [1 1; 1 lambda^(n-1)];
warning off;
m = m\[x(1); x(end)];
warning on;
a = m(1); b = m(2);

xi = a + b * lambda.^(0:n-1);
xi = min(xi, x(end));
xi = max(xi, x(1));

% interpolation
M1 = interp1( x, M, xi );
imagesc(x, y, M1);
set(gca, 'YScale', 'log');