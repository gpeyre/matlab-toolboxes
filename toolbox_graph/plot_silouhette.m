function plot_silouhette(M, s_list, light_pos, light_mode)

if nargin<4
    light_mode = 'local';
end
if nargin<3
    light_pos = [];
    light_mode = '';
end

% display
clf;
hold on;
C = ones(size(M(:,:,1)));
surf(M(:,:,1), M(:,:,2), M(:,:,3), C );
colormap gray(256);
for i=1:length(s_list)
    h = plot3(s_list{i}(1,:),s_list{i}(2,:),s_list{i}(3,:), 'r');
    set(h, 'LineWidth', 3);
end
% display the light position
if strcmp(light_mode, 'infinite')

else
    h = plot3(light_pos(1),light_pos(2),light_pos(3), 'b*');
    set(h, 'MarkerSize', 15);
end
hold off;

shading interp;
lighting phong;
camlight local; 
material shiny;
camproj('perspective');
axis square; axis off;
cameramenu;
view(3);
camlight;