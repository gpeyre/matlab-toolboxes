% test for optimization of the conditionning of the CS matrix

n = 100; m = 300;

[M,eOrth,eAngle,eNorm] = perform_cs_matrix_optimization(randn(n,m), options);


clf;
subplot(2,1,1);
plot(log(eOrth(3:end))); axis tight;
title('Orth');
subplot(2,1,2);
plot(log(eNorm(3:end))); axis tight;
title('Norm');

return;


%%% OLD CODE %%%

c = 0.15*2;
c = Inf;

niter = 20;
eOrth = [];
eNorm = [];
eAngle = [];
for i=1:niter
    % orthogonality constraint
    eOrth(end+1) = norm( M*M' - eye(n)*m/n );
    M = orth(M')';
    M = M(1:n,:) * sqrt(m/n);
    % norm constraint
    d = sqrt( sum( M.^2, 1 ) );
    eNorm(end+1) = sum( abs(d-1) );
    M = M./repmat( d, [n 1] );
    % angle constraint
    if 1
    G = M'*M;
    eAngle(end+1) = norm( G-eye(m) );
    G1 = min(max(G,-c),c);
    G1 = G1-diag(diag(G1))+diag(ones(1,m)); % remet ? 1 la diagonale
    [U,S,V] = svd(G1);
    M = U*diag(sqrt(diag(S)));
    M = U(:,1:n)';
    end
end


clf;
subplot(2,1,1);
plot(log(eNorm(3:end)));
subplot(2,1,2);
plot(log(eOrth(3:end)));