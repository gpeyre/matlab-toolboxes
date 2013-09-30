% test for div/grad operators

mydotp = @(x,y)sum(x(:).*y(:));

boundl = {'per', 'sym'};
for nbdims=2:3
    % randomized data
    if nbdims==2
        n = 128;
        f = randn(n);
        g = randn(n,n,2);
    else
        n = 32;
        f = randn(n,n,n);
        g = randn(n,n,n,3);
    end
    for i=1:2
        for order = 1:2
            options.order = order;
            options.bound = boundl{i};
            e = mydotp(grad(f,options), g) + mydotp(f, div(g,options));
            disp(['nb.dims=' num2str(nbdims) ', order=' num2str(order) ', bound=' options.bound ', error=(should be 0) ' num2str(e) '.']);
        end
    end
end