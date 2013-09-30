%% test for peano and hilbert curve


clf;
for i=1:4
    subplot(2,2,i);
    plot_curve( hilbert_curve(i), [], 'b' );
    title(sprintf('n=%d',i));
    axis off;
end
saveas(gcf,'courbe_hilbert','bmp');


clf;
for i=1:4
    subplot(2,2,i);
    plot_curve( peano_curve(i), [], 'b' );
    title(sprintf('n=%d',i));
    axis off;
end
saveas(gcf,'courbe_peano','bmp');