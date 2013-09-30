% test for plot of local DCT

n = 1024;

btree = [1   1 0   0 1 0 0   0 0 0 0 0 0 0 0];

options.btree = btree;

delta = 15;
a = [0 n/4+15 3/8*n n/2]+delta;


style = 'tree';
style = 'frame';
clf;
for i=1:length(a)
    subplot(length(a),1,i);
    x = perform_best_dct(dirac(n,a(i)),-1,options);
    plot_best_basis(x, btree, style);
end