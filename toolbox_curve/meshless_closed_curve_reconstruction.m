function I = meshless_closed_curve_reconstruction(c,R)

% meshless_closed_curve_reconstruction - reconstruct a closed curve from point cloud.
%
%   I = meshless_curve_reconstruction(c, R);
%
%   I is the ordering of the points in c (a 2 x npts matrix).
%   R is the locallity radius (the smaller, the sparser).
%
%   Copyright (c) 2004 Gabriel Peyré

npts = size(c,2);

x = c(1,:);
y = c(2,:);
if nargin<2
    R = abs(max(x)-min(y)) + abs(max(x)-min(y));
    R = R*20/npts;
end

[tmp,i] = min(x);

done = zeros(npts,1)+Inf;

last_ind = 0;

while ~isempty(done==Inf)
    
    % select a bunch of points arround i
    d = sqrt( ( x-x(i) ).^2 + ( y-y(i) ).^2 ); d=d';
    I = find( d<R & done==Inf);
    
    if isempty(I)
        break;  
    end
    if length(I)==1
        done(I) = last_ind+1;
        break;
    end
    
    if 0
        I1 = find(d<R & done==Inf);
        I2 = find(d>=R | done~=Inf);
        clf;
        hold on;
        plot(x(I1), y(I1), 'r.')
        plot(x(I2), y(I2), 'b.');
        hold off;
    end
    
    % project points on a plane
    xx = x(I);  yy = y(I);
    xs = xx-mean(xx);
    ys = yy-mean(yy);
    C = cov( [xs',ys'] );
    [U,S,V] = svd(C);
    v = V(:,1);
    s = xs*v(1)+ys*v(2);
    
    [tmp,Icur] = sort(s);
    
    done(I) = Icur + last_ind;
    last_ind = last_ind + length(Icur);
    
    i = I( Icur(end) );     % last point   
    
end

if 0
    I1 = find(done==Inf);
    I2 = find(done~=Inf);
    clf;
    hold on;
    plot(x(I1), y(I1), 'r.')
    plot(x(I2), y(I2), 'b.');
    hold off;
end
    
% in the other direction
[tmp,i] = min(done);
last_ind = 0;

while ~isempty(done==Inf)
    
    % select a bunch of points arround i
    d = sqrt( ( x-x(i) ).^2 + ( y-y(i) ).^2 ); d=d';
    I = find( d<R & done==Inf);
    
    if isempty(I)
        I = done;
        return; 
    end
    
    
    if 0
        I1 = find(d<R & done==Inf);
        I2 = find(d>=R | done~=Inf);
        clf;
        hold on;
        plot(x(I1), y(I1), 'r.')
        plot(x(I2), y(I2), 'b.');
        hold off;
    end
    
    % project points on a plane
    xx = x(I);  yy = y(I);
    xs = xx-mean(xx);
    ys = yy-mean(yy);
    C = cov( [xs',ys'] );
    [U,S,V] = svd(C);
    v = V(:,1);
    s = xs*v(1)+ys*v(2);
    
    [tmp,Icur] = sort(s);
    Icur = reverse(Icur);
    
    done(I) = -Icur + last_ind;
    last_ind = last_ind - length(Icur);
    
    i = I( Icur(end) );     % last point   
    
end


I = done;