function H = compute_symmetric_conditional_histogram(x, x_cond, p, pc)

% compute_symmetric_conditional_histogram - compute a symmetric histogram
%
%   H = compute_symmetric_conditional_histogram(x, x_cond, p, pc);
%
%   x is assumed to be integer valued.
%   p is the nbr of bins for value
%   pc is the nbr of bins for conditionning
%
%   Copyright (c) 2005 Gabriel Peyré

I = find( x>p | x<-p | x_cond>pc | x_cond<-pc );
x(I) = [];
x_cond(I) = [];

H = zeros(2*pc+1,2*p+1);

for i = -pc:pc
    I = find( x_cond==i );
    h = compute_histogram( x(I), -1, 0 );
    nh = (length(h)-1)/2;
    if nh>0
        H(i+pc+1, p+1-nh:p+1+nh) = h(:)';
    end
end
