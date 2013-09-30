function MW = perform_2d_int_translation(M, dx, dy, bound)

% perform_2d_int_translation - perform a 2D integer translation
%
%   MW = perform_2d_int_translation(M, dx, dy, bound);
%
%   bound tells how boundary are handeled, it can 
%   be either 'sym' (symmetric) or 'per' (periodic)
%
%   Copyright (c) 2005 Gabriel Peyré

if nargin<4
    bound = 'sym';
end

MW = M;
switch lower(bound)
    case 'sym'
        % translation along X
        if dx>0
            MW = [MW(1+dx:end,:); MW(end:-1:end-dx+1,:)];
        elseif dx<0
            MW = [MW(-dx:-1:1,:); MW(1:end+dx,:)];
        end
        % translation along Y
        if dy>0
            MW = [MW(:,1+dy:end), MW(:,end:-1:end-dy+1)];
        elseif dy<0
            MW = [MW(:,-dy:-1:1), MW(:,1:end+dy)];
        end
    case 'per'
        % translation along X
        if dx>0
            MW = [MW(1+dx:end,:); MW(1:dx,:)];
        elseif dx<0
            MW = [MW(end+dx+1:end,:); MW(1:end+dx,:)];
        end
        % translation along Y
        if dy>0
            MW = [MW(:,1+dy:end), MW(:,1:dy)];
        elseif dy<0
            MW = [MW(:,end+dy+1:end), MW(:,1:end+dy)];
        end
    otherwise
        error('Unknown boundary conditions');
end