% test for adaptive arithmetic coders
%
%   Copyright (c) 2006 Gabriel PeyrŽ

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bijectivity tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nb_coders = 4;
names = {'LetItWave', 'Escape', 'Skretting', 'Lepennec'};

% generate samples
n = 1000;
p = 100;
x = floor( rand(n,1)*p )-p/2;

for i=1:nb_coders
    disp(['--> Testing coder ' names{i}]);
    clear options;
    options.coder_type = i;
    %% storing data size in bit stream %%
    tic;
    y = perform_arithmetic_coding(x,1,options);
    xx = perform_arithmetic_coding(y,-1,options);
    toc;
    err = norm(x-xx, 'fro'); % should be 0
    disp( sprintf('     Error (should be 0): %.2f', err) );
    %% not storing the size of the data %%
    options.known_size = n;
    y = perform_arithmetic_coding(x,1, options);
    xx = perform_arithmetic_coding(y,-1, options);
    err = norm(x-xx, 'fro'); % should be 0
    disp( sprintf('     Error (should be 0): %.2f', err) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% efficiency tests

n_list = 100:200:7000;
% generate signal
options.evol = 0.5; % evolution in the probability
x = load_signal('rand', max(n_list), options);
h = compute_entropy(x);

nb_bits = [];
for i=1:nb_coders
    fprintf(['--> Testing coder ' names{i} ' ']);
    clear options;
    options.coder_type = i;
    nb = [];
    for n=n_list
        fprintf('.');
        [y,nb(end+1)] = perform_arithmetic_coding(x(1:n),+1,options);
    end
    fprintf('\n');
    nb_bits = [nb_bits nb(:)./n_list(:)];
end

clf;
plot(n_list, nb_bits-h, '.-');
axis tight;
disp('Bit per symbol - entropy');
disp('Size');
legend(names);