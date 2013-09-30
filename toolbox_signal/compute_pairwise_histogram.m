function H = compute_pairwise_histogram(f, x, xc, options)



if size(f,1)<size(f,2)
    f = f';
end
if size(x,1)<size(x,2)
    x = x';
end

nb_samples = Inf;
if isfield(options, 'nb_samples')
    nb_samples = options.nb_samples;
end
nb_samples = min(nb_samples,size(f,1));
sel = randperm(size(f,1)); sel = sel(1:nb_samples);
f = f(sel,:);

n = size(x,1);
nc = size(xc,1);

H = zeros(n,nc);

sc = xc(2)-xc(1);
s = x(2)-x(1);


for ic = 1:length(xc)
    I = find( abs(f(:,2)-xc(ic))<=sc/2 );
    fI = f(I,1);
    J = find( fI>min(x)-s/2 & fI<max(x)+s/2 );
    H(:,ic) = hist(fI(J),x);
end