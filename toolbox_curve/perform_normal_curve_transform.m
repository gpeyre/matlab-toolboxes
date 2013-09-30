function y = perform_normal_curve_transform(x,options)

% mc_forward - perform fwd multiresolution curve.
%
% Forward transform:
%   cw = perform_normal_curve_transform(c,options);
% Backward transform: 
%   c = perform_normal_curve_transform(cw);
%
%   options.type is the kind of predictor used, can be either 'linear' or 'cubic'
%   options.J is the number of levels of subdivision (ie. the transform will use 2^J coefs).
%
%   The wavelets coefficients are stored in 'cw.coef' and the extremity
%   points are in 'cw.begin' and 'cw.end'.
%
%   Copyright (c) 2006 Gabriel Peyré

options.null = 0;
if isfield(options, 'J');
    J = options.J;
else
    J = 10;
end
if isfield(options, 'type');
    type = options.type;
else
    type = 10;
end

if ~isstruct(x)
    %%% Forward %%%
    % c should be of size (2,n)
    if size(x,1)~=2
        x = x';
    end
    if size(x,1)~=2
        error('c is not a 2D curve.');
    end
    % take care of not having same end points
    if sum( abs(x(:,1)-x(:,end)) )<eps
        x(:,end) = [];
    end
    y = mc_forward( x, J, type );
else
    %%% Backward %%%
    y = mc_backward( x );
end


function c = mc_backward( mc )

% mc_backward - perform bwd multiresolution curve.
%
% c = mc_backward( mc )
%
%   Copyright (c) 2004 Gabriel Peyré


type = mc.type;

coef = mc.coef;
% nbr coef = 2^J-1;
J = floor( log2(length(coef)+1) );

% reconstructed curve at current level
c = [mc.begin, mc.end];


k = 1;  % number of extracted coef

for j=1:J
    
    P = size(c,2);  % number of points
    
    cc1 = [];   % new points to add at current level
    
    for i = 1:P-1
        
        x1 = c(:,i);
        x2 = c(:,i+1);
        
        % perform prediction
        if i>1
            x0 = c(:,i-1);
        else
            x0 = x1;
        end
        if i<P-1
            x3 = c(:,i+2);
        else
            x3 = x2;
        end
        
        % perform prediction
        if strcmp(type,'linear')
            x = (x1+x2)/2;    
        elseif strcmp(type,'cubic')
            x = -x0/16 + x1*9/16 + x2*9/16 - x3/16;
        else
            error('Unknown scheme.');
        end

        % compute normal
        n = x2-x1;  n = [-n(2); n(1)];  n = n / norm(n,'fro');
        
        y = x + n*coef(k);
        k = k+1;
        
        cc1 = [cc1,y];
        
    end
    
    % create new vector by entrelacing new/old
    cc_new = zeros(2,2*P-1);
    cc_new(:,1:2:end) = c;
    cc_new(:,2:2:end) = cc1;
    c = cc_new;    
end




function mc = mc_forward( c, J, type )

% mc_forward - perform fwd multiresolution curve.
%
% mc = mc_forward( c, J, type );
%
%   'type' is the kind of predictor used, can be either 'linear' or 'cubic'
%   J is the level of subdivision (ie. the transform will use 2^J coefs).
%
%   The wavelets coefficients are stored in 'mc.coef' and the extremity
%   points are in 'mc.begin' and 'mc.end'.
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<2
    J = 4;
end

if nargin<3
    type = 'linear';
end


mc.type = type;


mc.begin = c(:,1);
mc.end   = c(:,end);
mc.coef = [];

% reconstructed curve at current level
cc = [c(:,1), c(:,end)];
Icc = [1, size(c,2)];       % index of the selected points 

for j=1:J
    P = size(cc,2);  % number of points
    
    cc1 = [];   % new points to add at current level
    Icc1 = [];
    
    for i = 1:P-1
        
        x1 = cc(:,i);
        x2 = cc(:,i+1);
        Ix1 = Icc(i);
        Ix2 = Icc(i+1);
        
        if i>1
            x0 = cc(:,i-1);
        else
            x0 = x1;
        end
        if i<P-1
            x3 = cc(:,i+2);
        else
            x3 = x2;
        end
        
        % perform prediction
        if strcmp(type,'linear')
            x = (x1+x2)/2;    
        elseif strcmp(type,'cubic')
            x = -x0/16 + x1*9/16 + x2*9/16 - x3/16;    
        else
            error('Unknown scheme.');
        end

        % compute normal
        n = x2-x1;  n = [-n(2); n(1)];  
        if norm(n,'fro')>0
            n = n / norm(n,'fro');
        end
        % signed distance from point y in c to line (x,n) is d=(xy^u)/|u|
        xy = [ c(1,:) - x(1); c(2,:) - x(2) ];
        D = ( xy(1,:).*n(2)-xy(2,:).*n(1) );
        
        % locate zero crossing
        I1 = findzeros(D);
        if isempty(I1)
            warning('Fail to find an intersection point.');    
            [tmp, I1] = min(abs(D));
        end
        
        % extract the zeros that lies between x1 and x2
        if Ix1>Ix2 % swap index
            tmp = Ix1; Ix1 = Ix2; Ix2 = tmp;
        end
        
        I2 = find( I1>=Ix1 & I1<=Ix2 );
        if isempty(I2)
            warning('Fail to find an intersection point.');
            [m,I2] = min( abs(I1 - (Ix1+Ix2)/2) );
        end
        I = I1(I2(1));
        Icc1 = [Icc1, I];
        
        % find closest point
        y1 = c(:,I);  d1 = D(I);
        if I>1 && D(I)*D(I-1)<0
            y2 = c(:,I-1);    d2 = D(I-1);
        else
            y2 = c(:,I+1);    d2 = D(I+1);
        end
        % interpolate position
        y = abs(d2)/(abs(d1)+abs(d2))*y1 + abs(d1)/(abs(d1)+abs(d2))*y2;
        
        % the coefficient !
        mc.coef = [mc.coef, dot( y-x,n )];
        cc1 = [cc1,y];
        
    end
    
    % create new vector by entrelacing new/old
    cc(:,1:2:2*P-1) = cc;
    cc(:,2:2:2*P-1) = cc1;
    Icc(1:2:2*P-1) = Icc;
    Icc(2:2:2*P-1) = Icc1;
    
end




function varargout = findzeros(varargin)

%FINDZEROS  Find zeros in a vector.
%  IND=FINDZEROS(Y) finds the indices (IND) which are close to local zeros in the vector Y.  
%  [IND,YZEROS]=FINDZEROS(X,Y), besides IND finds the values of local zeros (YZEROS) in the 
%  vector X by linear interpolation of the vector Y.
%  Inputs:
%   X: vector
%   Y: vector
%  Outputs:
%   IND: indices of the zeros values in the vector Y
%   YZEROS: values in the vector X of zeros in the vector Y
%  Marcos Duarte  mduarte@usp.br 11oct1998 

if nargin==1
    y=varargin{1};
    if size(y,2)==1
        y=y';
    end
elseif nargin ==2
    x=varargin{1}; 
    y=varargin{2}; 
    if ~isequal(size(x),size(y))
        error('Vectors must have the same lengths')
        return
    elseif size(x,2)==1
        x=x'; y=y';
    end
else
    error('Incorrect number of inputs')
    return
end

% find +- transitions (approximate values for zeros):
%ind(i): the first number AFTER the zero value
ind = sort( find( ( [y 0]<0 & [0 y]>0 ) | ( [y 0]>0 & [0 y]<0 ) ) );
% find who is near zero:
ind1=find( ( abs(y(ind)) - abs(y(ind-1)) )<= 0 );
ind2=find( ( abs(y(ind)) - abs(y(ind-1)) )> 0 );
indzero = sort([ind(ind1) ind(ind2)-1 find(y==0)]);
varargout{1}=indzero;
if exist('x')
    % find better approximation of zeros values using linear interpolation:
    yzeros=zeros(1,length(ind));
    for i=1:length(ind)
        p=polyfit( y(ind(i)-1:ind(i)),x(ind(i)-1:ind(i)),1 );
        yzeros(i)=polyval(p,0);
    end
    yzeros=sort([yzeros x(find(y==0))]);
    varargout{2}=yzeros;
end