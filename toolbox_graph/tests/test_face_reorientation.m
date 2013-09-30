% test for mesh re-orientation

path(path, '../toolbox_graph/off');
rep = '../toolbox_graph/off/';
name = 'david50kf';
name = 'nefertiti';
name = 'beetle';
name = 'fandisk';
[vertex,faces] = read_mesh(name);

% flip faces at random
n = size(vertex,2);
m = size(vertex,2);
I = find(rand(m,1)>.5);
faces(1:2,I) = faces(2:-1:1,I);


options.method = 'fast';
options.method = 'slow';
faces1 = perform_faces_reorientation(vertex,faces,options);


clf;
plot_mesh(vertex,faces1);
shading flat;