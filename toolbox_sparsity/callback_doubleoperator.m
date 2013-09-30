function y = callback_doubleoperator(x,dir,options)

% callback_doubleoperator - apply two calbacks
%
%   y = callback_doubleoperator(x,dir,options);
%
%   Copyright (c) 2008 Gabriel Peyre

c1 = getoptions(options, 'callback1', 0, 1);
c2 = getoptions(options, 'callback2', 0, 1);



if dir==1
    y = feval(c1, x, dir, options);
    y = feval(c2, y, dir, options);
else    
    y = feval(c2, x, dir, options);
    y = feval(c1, y, dir, options);
end