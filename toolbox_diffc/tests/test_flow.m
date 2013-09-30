% test for flow creation
% create a flow with a vortice and a source.

n = 70;
p = [[0.2;0.2], [0.8;0.8]];
r = [1,0];  % rotation
s = [0,1];  % source
c = [2,0];

V = compute_flow(p, r, s, c, n);

Vn = perform_vf_normalization(V);
Vc = clamp(V, -10,10);  % clamp for display
clf;
plot_vf(Vn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute streamlines
x = 0:1/(n-1):1;
[X,Y] = meshgrid(x,x);

p = 20;
sx = 0:1/(p-1):1;
sy = 0*sx + 1/2;
% up
XY = stream2(X,Y,Vc(:,:,2),Vc(:,:,1),sx,sy);
% down
XY1 = stream2(X,Y,-Vc(:,:,2),-Vc(:,:,1),sx,sy);
% merge
XY = {XY{:},XY1{:}};
% reverse stream
for i=1:length(XY)
    XY{i} = XY{i}(:,2:-1:1);
end

clf;
hold on;
plot_vf(Vn);
streamline(XY');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use matlb builting stream placement
[XY,tmp] = streamslice(X,Y,Vc(:,:,2),Vc(:,:,1), 2);
% reverse stream
for i=1:length(XY)
    XY{i} = XY{i}(:,2:-1:1);
end


clf;
hold on;
plot_vf(Vn, [], struct('linestyle', 'r'));
streamline(XY);
hold off;

