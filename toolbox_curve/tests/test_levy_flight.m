% test for 2D levy flight generation

save_image = 0;
n = 1024*4;
% exponent
alpha = 1;
% variance
sigma = 1;
% mean/median
beta = 0;
delta = 0;

type = 'axis';
type = 'isotropic';

alpha_list = linspace(0.6, 2, 9);
m = length(alpha_list);

clf;
for i=1:m
    alpha = alpha_list(i);
    x = gen_levy_flight(n,alpha,sigma,beta,delta,type);
    subplot(sqrt(m),sqrt(m), i);
    plot_curve(x);
    axis off;
    title( sprintf('\\alpha=%.1f', alpha) );
end

if save_image
    saveas(gcf, 'levy_flights.png');
end