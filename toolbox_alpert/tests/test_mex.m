% test for mex file speed
%
%   Copyright (c) 2004 Gabriel Peyré

% clear;
test_reconstruction = 0;
test_moments = 0;
test_speed = 1;
test_accuracy = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test for construction / reconstruction
if test_reconstruction
    n = 100;
    pos = rand(2,n);
    v = rand(n,1);
    monomials = [ [0;0], [0;1], [1;0], [2;0], [0;2], [1;1] ];
    % monomials = [ [0;0], [1;0] ];
    
    w = perform_moment_transform(v,pos, monomials, 1);
    v1 = perform_moment_transform(w,pos, monomials, -1);
    norme(v-v1)     % should be 0
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test for various moment configuration
if test_moments
    n = 150;
    for k=1:4
        
        pos = rand(2,n);
        v = rand(n,1);
        
        options.use_mex = 1;
        [w_mex,info] = perform_alpert_transform_1p5d(v,pos,k, 1, options);
        options.use_mex = 0;
        [w_mat,info] = perform_alpert_transform_1p5d(v,pos,k, 1, options);
        
        sw_mex = rev_sort_abs(w_mex);
        sw_mat = rev_sort_abs(w_mat);
        plot(1:n, sw_mex, 'r.', 1:n, sw_mat);
        legend('mex', 'mat');
        n_mex = norme(sw_mex);
        n_mat = norme(sw_mat);
        n_ori = norme(v);
        disp( sprintf('diff = %.3f, cons.mex=%.3f, cons.mat=%.3f', norme(sw_mex-sw_mat), abs(n_ori-n_mex), abs(n_ori-n_mat) ) );
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test for speed
if test_speed
    n = 1000;
    nb_test = 100;
    v = rand(n,1);
    pos = rand(2,n);
    
    monomials = [ [0;0], [0;1], [1;0], [2;0], [0;2], [1;1] ];
%    monomials = [ [0;0], [1;0] ];
    k = monomials;
    clear options;
    options.degree_type = 'user_defined';
    
    options.use_mex = 1;
    tic;
    for i = 1:nb_test
        [w_mex,info] = perform_alpert_transform_2d(v,pos,k, 1, options);
    end
    disp( ['time mex : ', num2str(toc), 's' ] );
    
    options.use_mex = 0;
    tic;
    for i = 1:nb_test
        [w_mex,info] = perform_alpert_transform_2d(v,pos,k, 1, options);
    end
    disp( ['time non-mex : ', num2str(toc), 's' ] );
    
    tic;
    for i = 1:nb_test
        w = perform_moment_transform(v,pos, monomials, 1);
    end
    disp( ['time mex+partition : ', num2str(toc), 's' ] );
    
    options.use_mex = 0;
    tic;
    for i = 1:nb_test
        w_haar = perform_haar_transform(v);
    end
    disp( ['time haar : ', num2str(toc), 's' ] );
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test for accuracy
if test_accuracy
    n = 1000;
    pos = (2*rand(2,n)-1);
    x = pos(1,:);
    y = pos(2,:);
    v = cos(x.^2 + y.^2)';

    sigma= 0.1;
    v = x.*exp(-x.^2/sigma); % +sin(y/5);
    
    monomials = [ [0;0], [0;1], [1;0], [2;0], [0;2], [1;1] ];
    monomials = [ [0;0], [1;0] ];
       
    % pure mex
    clear options;
    options.ptmax = size(monomials,2)/2;
    w1 = perform_moment_transform(v,pos, monomials, 1);
    
    % mex except for part
    options.use_mex = 1;
    k = monomials;
    options.degree_type = 'user_defined';
    [w2,info] = perform_alpert_transform_2d(v,pos,k, 1, options);
    
    % haar transform
    [tmp,I] = sort(pos(1,:));
    w3 = perform_haar_transform(v(I));
    
    
    loglog(1:n, rev_sort_abs(w1), 1:n, rev_sort_abs(w2), '.');
    loglog(1:n, l2error(w1), 1:n, l2error(w2), 'r--', 1:n, l2error(w3), 'g--');
    legend('mex-part', 'mex', 'haar');
end