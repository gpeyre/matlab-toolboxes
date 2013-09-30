function v = compute_rbf(x,d,xi, name)

switch name
    case 'abs'
        f = @(x)abs(x);
    case 'gauss'
        sigma = .03;
        f = @(x)exp(-x.^2/(2*sigma^2));
    case 'poly3'
        q = 3;
        f = @(x)abs(x).^q;
    case 'sqrt'
        q = .5;
        f = @(x)abs(x).^q;
    case 'thinplate'
        f = @(x)(x.^2).*log(abs(x)+.001);
end

n = length(xi);
m = length(x);

[Y,X] = meshgrid(x,x);

D = f(X-Y);

% solve for weights
a = pinv(D)*d;

v = f( repmat(xi, [1 m]) - repmat(x', [n 1]) ) .* repmat( a', [n 1] );
v = sum(v, 2);
