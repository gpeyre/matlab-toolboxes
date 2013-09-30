% test for hufman coding

% probability table
p = [.1 .15 .4 .15 .2];
m = length(p);

m = 16;
p = rand(m,1);
p = p/sum(p);

T = compute_hufftree(p);

% randomized signal
n = 1024;
x = round( rand(n,1)*(m-1) ) + 1;

y = perform_huffcoding(x,T,+1);

xx = perform_huffcoding(y,T,-1);

err = norm(x-xx);
disp(['Error (should be 0)=' num2str(err)'.']);

%%
% Display the computed tree
clf;
plot_hufftree(T);