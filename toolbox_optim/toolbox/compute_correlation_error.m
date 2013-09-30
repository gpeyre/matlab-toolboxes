function [err,ratio] = compute_correlation_error(Historic,fSolution, Gmat,F,q)

H = Historic./repmat(sqrt(sum(Historic.^2)), [length(fSolution) 1]);
S = repmat(fSolution, [1 size(Historic, 2)])/norm(fSolution);
err = 1 - sum(S.*H);
err = err(:);

ratio = [];
if nargin>2
    p = size(Gmat,1)/q;
    G  = @(f)reshape(Gmat*f, [p q]);
    Amplitude = @(u)sqrt(sum(u.^2,2));
    E = @(u)sum(Amplitude(u),1);
    J = @(f)E(G(f));
    R = @(f)-sum(F(:).*f) / J(f);
    r0 = R(fSolution);
    for i=1:size(Historic,2)
        progressbar(i,size(Historic,2));
        ratio(end+1) = R(Historic(:,i)) - r0;
    end
    ratio = ratio(:);
end