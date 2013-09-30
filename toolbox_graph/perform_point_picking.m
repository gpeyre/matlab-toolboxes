function  perform_point_picking( PointCloud, face )

% function    perform_point_picking( PointCloud, face );
%
%   This function shows a 3D point cloud or a mesh and lets the user click select one 
%   of the points by clicking on it. The selected point will be highlighted 
%   and its index in the point cloud will be printed. 
%
% input:
%   PointCloud:  should be a 3*N matrix representing N 3D points
%                (if it's M*N, M > 3, the higher dimensions are ignored)
%
% output:
%   none
%       
%   To test this function, run these lines:
%
%                 PointCloud = rand(3,100)*100;
%                 perform_point_picking( PointCloud );
%
%   now rotate or move the point cloud and try it again.
%   (on the figure View menu, turn the Camera Toolbar on, ...)
%
% by Babak Taati
% http://qlink.queensu.ca/~3bt1/
% Queen's University
% May 4, 2005 (revised Oct 30, 2007)


if nargin==2
    clf;
    hold on;
    plot_mesh(PointCloud,face);
end
plot3(PointCloud(1,:), PointCloud(2,:), PointCloud(3,:), 'c.'); % visualize the point cloud
hold on; % so we can highlight the clicked points without clearing the point cloud
cameramenu;

global PointCloudIndex;
set(gcf,'WindowButtonDownFcn',{@callbackClickA3DPoint,PointCloud}); % set the callbak
set(gcf,'WindowButtonUpFcn',{}); % set the callbak


% callbackClickA3DPoint.m
% by Babak Taati
% http://qlink.queensu.ca/~3bt1/
% Queen's University
% May 4, 2005
%
%
% this is the callback function for ClickA3DPoint
%   

function callbackClickA3DPoint(src, eventdata, PointCloud);

global PointCloudIndex;
Point = get(gca,'CurrentPoint');        % get the mouse click position

CamPos = get(gca,'CameraPosition');     % camera position
CamTgt = get(gca,'CameraTarget');       % where the camera is pointing to

CamDir = CamPos - CamTgt;               % Camera direction

CamUpVect = get(gca,'CameraUpVector');  % camera 'up' vector


% build a orthonormal frame based on the viewing direction and the up vector (the "view frame")
Z_axis = CamDir/norm(CamDir);    
UP_axis = CamUpVect/norm(CamUpVect); 
X_axis = cross(UP_axis,Z_axis);
Y_axis = cross(Z_axis, X_axis);


Rot = [ X_axis ; Y_axis ; Z_axis ];     % view rotation 

RotatedPointCloud = Rot * PointCloud;   % the point cloud represented in the view frame
RotatedPointFront = Rot * Point' ;      % the clicked point represented in the view frame


% --- find the nearest neighbour to the clicked point 
%
% project the 3D point cloud onto a plane (i.e. ignore the z coordinate) 
% and find the point that is the closest to the clicked point in that 2D plane.


diff = [ RotatedPointCloud(1,:) - RotatedPointFront(1)   ;  RotatedPointCloud(2,:) - RotatedPointFront(2) ]; % the difference of all the 2D points from the clicked point (2D)
AllDistances = RowNorm(diff');                  % the distance of all the 2D points from the clicked point (2D)
[dist, PointCloudIndex] = min(AllDistances);	% find the point that is closest to the clicked point (2D)


% ---- delete the old nearest neighbour & display the new nearest neighbour

h   = findobj(gca,'Tag','pt'); % try to find the old point
SelectedPoint = PointCloud(:,PointCloudIndex); 

if isempty(h) % Check if it's the first click (i.e. no need to delete the previous point)
    
    h = plot3(SelectedPoint(1,:), SelectedPoint(2,:), SelectedPoint(3,:), 'r.', 'MarkerSize', 20); % highlight the selected point
    set(h,'Tag','pt');   % set its Tag property for later use   

else    % if Not the first click

    delete(h);      % delete the previously selected point
    
    h = plot3(SelectedPoint(1,:), SelectedPoint(2,:), SelectedPoint(3,:), 'r.', 'MarkerSize', 20);  % highlight the newly selected point
    set(h,'Tag','pt');  % set its Tag property for later use


end


fprintf('you clicked on point number %d\n', PointCloudIndex);


% RowNorm.m
% Babak taati
%
% function  N = RowNorm(A);
%
% returns a matrix with same number of rows as A and one column
% each element in N is the norm of the corresponding row in A
%
% input:    A (an m*n matrix)
% output:   N (an m*1 vector)
%

function N = RowNorm(A);

N = sqrt(sum(A.*A, 2));

