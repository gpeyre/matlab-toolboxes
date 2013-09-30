% test the spherical compression scheme

sampling_type = 'random';
% sampling_type = 'uniform';

data_type = 'syntetic';
% data_type = 'real';

n = 1500;
k = 2;          % number of vanishing moments.
dir = 1;

switch sampling_type
    case 'random'
        pos = randn(3,n);
        d = sqrt( pos(1,:).^2+pos(2,:).^2+pos(3,:).^2 );
        pos(1,:) = pos(1,:)./d;
        pos(2,:) = pos(2,:)./d;
        pos(3,:) = pos(3,:)./d;
    case 'uniform'
        [pos,L,diam] =  sphere_sampling(n);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute some function
switch data_type
    case 'syntetic'
        v = pos(1,:).^2 - pos(2,:).*pos(3,:);
        v = v(:);
        v = rescale(v);
        v = cos(10*v);
        % list for compression plot
        m_list = [1,2,5,10,15];
    case 'real'
        sampling_type = 'real';
        load('data/data_lenglet.mat');
%        [pos,face] = gen_base_mesh('oct');
%        [pos,face] = subdivide_sphere(pos,face,2); pos = pos'; 
        v = ADC24_36_14(:);
        pos = grad_81';
        % symmetrize the data
        pos = [pos,-pos];
        v = [v;v];
        n = length(v);
        v = rescale(v);
        % list for compression plot
        m_list = [5,10,20,40,70];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute clustering
clear options;
options.ptmax = k^2;
options.part_type = 'kmeans';
[part,B,G] = dichotomic_grouping(pos,options);
% plot
clf;
subplot(1,2,1);
plot_spherical_partition(pos,part);
title('Clustering of sampling locations');

subplot(1,2,2);
plot_spherical_function(pos,v);
title('Function to transform (interpolated)');
saveas(gcf, ['function_',data_type, '_clustering_', sampling_type], 'png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform transform
clear options;
options.part_type = 'kmeans';
[w,info,part] = perform_alpert_transform_nd(v,pos,k, 1, options);
% try reconstruction just to check the transform is bijective ...
options.part = part;
v1 = perform_alpert_transform_nd(w,pos,k, -1,options);
norm(v-v1, 'fro') % should be 0

cv = rev_sort_abs(v);
cw = rev_sort_abs(w);
cv = l2error(v);
cw = l2error(w);

a = (1:n/2)';
x = log2(a);
y1 = cv(a);
y2 = cw(a);

% reglin
fiton = a;
A = a;
B = y2;
aa = log2(A(fiton));
bb = log2(B(fiton));
[coef,S] = polyfit(aa,bb,1);
disp( sprintf('----> 1D Alpert transform, decreasing=%.4f', coef(1)) );
reglin = log2(a)*coef(1)+coef(2);


clf;
hold on;
plot( x, log2(y1), x, log2(y2), x, reglin, 'b:' );
axis tight;
legend('Original error decreasing', 'Transformed error decreasing');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform reconstruction
clf;
subplot(2,3,1);
plot_spherical_function(pos,v);
title('Original function');
s = 1;

for m = m_list
    s = s+1;
    w1 = keep_biggest(w, round(n*m/100) );
    % backward transform
    v1 = perform_alpert_transform_nd(w1,pos,k, -1,part);
    subplot(2,3,s);
    plot_spherical_function(pos,v1);
    err = round( 100*norm(v-v1, 'fro')^2/norm(v, 'fro')^2 );    % RMS error in %
    title(sprintf('%.1f %% of coef, err=%d%%', m, err  ));
end

saveas(gcf, ['compression_',data_type, '_clustering_', sampling_type], 'png');