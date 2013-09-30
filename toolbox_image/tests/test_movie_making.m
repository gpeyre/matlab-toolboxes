% test for movie making

n = 128;
m = 128;
x = linspace(0,1,n^2);
M = zeros(n,n,m);
M(:,:,1) = perform_histogram_equalization(randn(n),x);
for i=2:m
    progressbar(i,m);
    M(:,:,i) = perform_blurring(M(:,:,i-1), 5);
    M(:,:,i) = perform_histogram_equalization(M(:,:,i), x);
end

rep = 'results/movies/';
if not(exist(rep))
    mkdir(rep);
end

compute_movie_file(M, [rep 'test.avi']);
compute_movie_file(M, [rep 'test.gif']);