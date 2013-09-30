% test for the haar transform

n = 64;
m = 1000;
J = 5;
x = rand(n,m);

y = perform_haar_transform(x,+1,J);
x1 = perform_haar_transform(y,-1,J);

% check reconstruction
err = sum( (x(:)-x1(:)).^2 );
disp(['Reconstruction error (shoule be 0) : ' num2str(err, 3)]);
% check energy conservation
E1 = sum( x.^2, 1 );
E2 = sum( y.^2, 1 );
err = sum( abs(E1(:)-E2(:)) );
disp(['Energy conservation (shoule be 0) : ' num2str(err/sum(abs(E1(:)))*100, 3) '%']);


% test for translation invariant
n = 512;
options.ti = 1;
x = load_signal('Piece-Regular', n);
J = 5;
y = perform_haar_transform(x,+1,J, options);
x1 = perform_haar_transform(y,-1,J, options);

d = reshape(y,n,J+1);

E1 = sum( x.^2, 1 );
E2 = sum( sum( d.^2 ) ./ 2.^[J J:-1:1] );
err = sum( abs(E1(:)-E2(:)) );
disp(['Energy conservation (shoule be 0) : ' num2str(err/E1*100, 3) '%'] );

err = sum( (x(:)-x1(:)).^2 );
disp(['Reconstruction error (shoule be 0) : ' num2str(err, 3)]);