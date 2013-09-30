% test for the pyramidal decomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.bound = 'per';
% options.bound = 'sym';
options.J = 5;

g = MakeONFilter('Haar',1);
%g = MakeONFilter('Battle',3);

[qmf,dqmf] = MakeBSFilter( 'CDF', [4,4] );
g = {qmf,dqmf};

sigma = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test in 1D
disp('--- test in 1D ---');
n = 1024;
x = (0:1/(n-1):1)';
% x = ones(n,1);
y = cos( 8*pi*x.^2 );
y = x.^2;

y = y + randn(size(y))*sigma;

% fwd transform
f = perform_pyramid_transform(y,g,options);

% threshold
N = 100;
f = keep_biggest(f,N);

% bwd transform
y_fr = perform_pyramid_transform(f,g,options);
y_nfr = perform_pyramid_transform_nonframe(f,g,options);

Eo = norme(y)^2;
Ep = 0;
for i=1:length(f)
    Ep = Ep + norme(f{i})^2;
end

disp(sprintf('Energy : original=%.3f, transformed=%.3f.', Eo, Ep));

disp( sprintf('Error frame=%.3f, Error non-frame=%.3f', norme(y-y_fr), norme(y-y_nfr) ) ); 
disp( sprintf('PSNR frame=%.3f, PSNR non-frame=%.3f', psnr(y,y_fr), psnr(y,y_nfr)) ); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test in 2D
disp('--- test in 2D ---');
n = 512;
y = load_image('peppers', n);
options.J = 6;

% fwd transform
f = perform_pyramid_transform(y,g,options);

% threshold
if 1
N = 2^14;
f = keep_biggest(f,N);
end

% bwd transform
y_fr = perform_pyramid_transform(f,g,options);
y_nfr = perform_pyramid_transform_nonframe(f,g,options);

disp( sprintf('Error frame=%.3f, Error non-frame=%.3f', mean(mean((y-y_fr).^2)), mean(mean((y-y_nfr).^2)) ) ); 
disp( sprintf('PSNR frame=%.3f, PSNR non-frame=%.3f', psnr(y,y_fr), psnr(y,y_nfr)) ); 
