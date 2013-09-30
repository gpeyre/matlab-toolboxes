% plot of sampling location on warped grid


rep = 'results/';

if ~exist(rep)
    mkdir(rep);
end

save_warped = 0;
save_dich = 1;

s = 50;
c = 0.45;
n = 16;
t = linspace(0,1,n);
[Y,X] = meshgrid(t,t);

x = X(:); y = Y(:);

xw = x; 
yw = y + c * x.^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot of the warped locations
clf;
hold on;
scatter(x,y, s, 'filled');
plot(t, 1-c*t.^2);
hold off;

clf;
scatter(x,y, s, 'filled');
axis equal;
axis off;
if save_warped
saveas(gcf, [rep 'samples.eps'], 'epsc');
end

clf;
scatter(xw,yw, 25, 'filled');
axis equal; axis off;
if save_warped
saveas(gcf, [rep 'samples_warped.eps'], 'epsc');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot of the grouping process
posw = [xw, yw]';
pos = [x, y]';
options.part_type = '1axis';
for i=0:5
    options.depthmax = i;
    [part,B,G] = dichotomic_grouping([yw, xw]',options);
    options.point_size = 80;
    clf;
    plot_dichotomic_partition(posw, part, options);
    axis([0 1.5 0 1.5]); axis equal; axis off;
    if save_dich
        saveas(gcf, [rep '_dich_' num2str(i) '.eps'], 'epsc');
    end
    clf;
    plot_dichotomic_partition(pos, part, options);
    axis([0 1.5 0 1.5]); axis equal; axis off;
    if save_dich
        saveas(gcf, [rep '_dich_' num2str(i) 'w.eps'], 'epsc');
    end 
end