% test for normal multiresolution curve

n = 600;
name = 'cloudy';
name = 'spiral';
name = 'levy';
name = 'rand';
options.alpha = 0.5;
c = load_curve(name,n,options);

options.J = 10;
options.type = 'linear';
mcL = perform_normal_curve_transform( c, options );
options.type = 'cubic';
mcQ = perform_normal_curve_transform( c, options );

clf;

i = 0;
for k = [10,15,25,100]
    
    i = i+1;
    mc1 = mcL;
    mc1.coef = keep_biggest(mc1.coef,k);
    subplot(2,4,i);
    cc = perform_normal_curve_transform( mc1 );
    plot_curve( cc );
    axis tight;
    axis equal;
    axis off;
    title(sprintf('%d coefs', k));

    mc1 = mcQ;
    mc1.coef = keep_biggest(mc1.coef,k);
    subplot(2,4,i+4);
    plot_curve( perform_normal_curve_transform( mc1 ) );
    axis tight;
    axis equal;
    axis off;
end
pause;

% decreasing of the coefficients
P = length(mcL.coef);
clf;
loglog( 1:P, rev_sort_abs(mcL.coef), 1:P, rev_sort_abs(mcQ.coef)  );
legend('linear', 'cubic');
axis tight;