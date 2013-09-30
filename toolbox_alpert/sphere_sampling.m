function [Points,L,diam,topcol,botcol] = sphere_sampling(N,lrounded,angles,mrounded,ipl)

% sphere_sampling - sample points on a sphere.
%
% [Points,L,diam] =  sphere_sampling(N,lrounded,angles,mrounded,ipl)

%   $Revision: 1.1 $ Paul Leopardi 2003-10-13
%   Make angles=1 the default
%   $Revision: 1.1 $ Paul Leopardi 2003-09-18
%   Correct calculation of diameters
%   Return topcol and botcol: number of regions in top and bottom collars
%   $Revision: 1.1 $ Paul Leopardi 2003-09-13
%   Parameter angles controls whether to locate points via angle or height average
%   Match twist and offset to Maple version of the algorithm
%   $Revision: 1.1 $ Paul Leopardi 2003-09-09
%   Parameters lrounded and mrounded control rounding of L and m
%   $Revision: 1.1 $ Paul Leopardi 2003-09-08
%   Return diameters
%   $Revision: 1.1 $ Paul Leopardi 2003-09-08
%   for UNSW School of Mathematics

global area;
if nargin < 2
    lrounded = 0;
end
if nargin < 3
    angles = 1;
end
if nargin < 4
    mrounded = 0;
end
if nargin < 5
    ipl = 0;
end
if N == 1
    Points = zeros(3,N);
    Points(:,1)=[0,0,1]';
    L = 1;
    diam = 2*pi;
    topcol = 0;
    botcol = 0;
else
    area = 4*pi/N;
    % [A,val,exitflag,output] = fsolve(@square,sqrt(area),optimset(optimset('fsolve'),'Display','off'));
    % A = real(A);
    A = sqrt(area);
    Beta = acos(1-2/N);
    fuzz = eps*2*N;
    if lrounded == 1
        L = 2+max(1,round((pi-2*Beta)/A));
    else
        L = 2+max(1,ceil((pi-2*Beta)/A-fuzz));
    end
    Theta = (pi-2*Beta)/(L-2);
    mbar = zeros(1,L);
    mbar(1) = 1;
    for i = 2:L-1
        mbar(i) = N*(cos(Theta*(i-2)+Beta)-cos(Theta*(i-1)+Beta))/2;
    end
    mbar(L) = 1;
    alpha = 0;
    m = zeros(1,L);
    m(1) = 1;
    diam = zeros(1,L);
    diam(1) = 2*Beta;
    for i = 2:L
        if mrounded == 1
            m(i) = round(mbar(i)+alpha);
        else
            if (mbar(i)-floor(mbar(i)+alpha+fuzz)) < 0.5
                m(i) = floor(mbar(i)+alpha+fuzz);
            else
                m(i) = ceil(mbar(i)+alpha-fuzz);
            end
        end
        alpha = alpha + mbar(i)-m(i);
    end
    
    if ipl > 0
        fprintf('For N=%d, the number of levels cuts L=%d\n', N, L);
        fprintf('--------------------------------------------------\n');
    end
    Points = zeros(3,N);
    Points(:,1) = [0,0,1]';
    z = 1-(2+m(2))/N;
    Format = 'Level %3d: mbar(%2d) =%12.8f   m(%2d) =%4d diff=%12.8f; Area= %9.5f\n';
    i = 1;
    Area = 1;
    if ipl > 0
        fprintf(Format,i,i,mbar(i),i,m(i),mbar(i)-m(i),Area);
    end
    offset = zeros(1,L-1);
    offset(1) = 0;
    n=1;
    for i = 2:L-1
        twist=4;
        if m(i-1) ~= 0 & m(i) ~= 0
            offset(i) = offset(i-1)+gcd(m(i),m(i-1))/(2*m(i)*m(i-1))+min(twist,floor(m(i-1)/twist))/m(i-1);
        else
            offset(i) = 0;
        end
        top = z+m(i)/N;
        rtop = sqrt(1-top^2);
        bot = z-m(i)/N;
        rbot = sqrt(1-bot^2);
        if m(i) > 1
            angle = 2*pi/m(i);
            if rtop > rbot
                rside = rtop;
                hside = top;
            else
                rside = rbot;
                hside = bot;
            end
            side = acos([rside,0,hside]*[rside*cos(angle),rside*sin(angle),hside]');
            diag = acos([rtop,0,top]*[rbot*cos(angle),rbot*sin(angle),bot]');
            diam(i) = max(side,diag);
        else
            diam(i) = acos([rtop,0,top]*[-rbot,0,bot]');
        end
        if angles
            h = cos((acos(top)+acos(bot))/2);
        else
            h = z;
        end
        r = sqrt(1-h^2);
        for j = 0:m(i)-1
            s = offset(i)+j/m(i);
            n = n+1;
            Points(:,n) = [r*cos(2*pi*s),r*sin(2*pi*s),h]';
        end
        Area = 1;
        if ipl > 0
            fprintf(Format,i,i,mbar(i),i,m(i),mbar(i)-m(i),Area);
        end
        z = z-(m(i)+m(i+1))/N;
    end
    i = L;
    Area = 1;
    if ipl > 0
        fprintf(Format,i,i,mbar(i),i,m(i),mbar(i)-m(i),Area);
    end
    Points(:,N) = [0,0,-1]';
    diam(L) = 2*Beta;
    if L < 3
        topcol = 0;
        botcol = 0;
    else
        topcol = m(2);
        botcol = m(L-1);
    end
end

function v=square(A)
global area;
v=8*asin(sqrt(2)*sin(A/2)/sin(A))-2*pi-area;
