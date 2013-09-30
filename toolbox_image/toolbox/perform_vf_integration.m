function M = perform_vf_integration(vf, dt,T_list, options )

% perform_vf_integration - perform a time integration of the vf
%
%   M = perform_vf_integration(vf, dt, T_list, options );
%
%   If you set options.posx or options.posy, then 
%       M(i,j,:,t)) is the (X,Y) location of the point (options.posx(i),options.posy(j))
%       at timestep T_list(t).
%   If you set options.pos, then M(i,:,t) is the location of the point
%       options.pos(:,i).
%
%   The vf is assumed to be sampled at locations [1,...,n]^2 unless posx and posy are given.
%
%   The integrator is a simple explicit euler, but it is very 
%   fast thanks to vectorization.
%
%   Copyright (c) 2005 Gabriel Peyre

n = size(vf,1);
p = size(vf,2);
m = length(T_list);

options.null = 0;

bound = getoptions(options, 'bound', 'sym');
verb = getoptions(options, 'verb', 0);

if isfield(options, 'pos')
    pos = options.pos;
    posx = [];
else
    posx = getoptions(options, 'posx', 1:n);
    posy = getoptions(options, 'posy', 1:p);
    [Y,X] = meshgrid(posy,posx);
    pos = cat(1,X(:)',Y(:)');
end

% number of points
q = size(pos,2);

% M = zeros( length(posx),length(posy),m);
M = zeros(q,2,m);

T_list = [0; sort(T_list(:))];

for i=1:m
    if verb
        progressbar(i,m);
    end
    % compute the number of time steps
    delta = T_list(i+1)-T_list(i);
    nt = ceil(delta/dt);
    if nt>0
        DT = delta / nt;
        % perform each time step
        for k=1:nt
            dx = interp2(1:p,1:n,vf(:,:,1), pos(2,:),pos(1,:) );
            dy = interp2(1:p,1:n,vf(:,:,2), pos(2,:),pos(1,:) );
            pos = pos + DT*cat(1,dx,dy);
%            X = X + DT*dx; Y = Y + DT*dy;
            switch bound
                case 'sym'
                    pos(1,:) = clamp( pos(1,:), 1,n );
                    pos(2,:) = clamp( pos(2,:), 1,p );
                case 'per'
                    pos(1,:) = mod(pos(1,:)-1,n-1)+1;
                    pos(2,:) = mod(pos(2,:)-1,p-1)+1;
                otherwise
                    error('Unknown options.bound.');
            end
        end
    end
    % record value
    M(:,:,i) = pos'; % X + 1i * Y;
    % check for early break if nothing is moving anymore    
    if i>1
        d = M(:,:,i) -  M(:,:,i-1);
        d = mean( sqrt(sum(d.^2,2)) );
        if d<1e-6
            progressbar(m,m);
            break;
        end
    end
end
M = M(:,:,1:i);
if not(isempty(posx))
    M = reshape(M, [length(posx) length(posy) 2 m]);
end
