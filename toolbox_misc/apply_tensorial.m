function y = apply_tensorial(x,func_name,func_arg)

% apply_tensorial - Multimensional application of a 1D function.
%
%   Apply a given 1D function to each dimension of the array x.
%   The name of the function is given in the string 'func_name'.
%   This is very powerful because it always send *row* vectors to the 
%   function, never vector with size e.g. [1,1,4,1,1] ...
% 
%   You can give additional arguments in the string func_arg.
%   Suppose you want to apply function my1Dfunc(x,'This is a string argument') to 
%   each dimension of a matrix X :
%       X = rand(3,3,3);
%       Y = apply_tensorial(X, 'my1Dfunc', '''This is a string argument''');
%
%   Copyright (c) 2003 Gabriel Peyré

y = x;
dim = size(x);
d = length(dim);

if( nargin==2 )
    func_arg='';
end

for s=1:d
    % perform a 1D transform in the dimension s
    % the remainging dimensions
    curd = dim([1:(s-1),s+1:d]);
    for t=1:prod(curd)
        % turn t into an index
        [out{1:d-1}] = ind2sub(curd,t);
        ind = cell2mat(out);
        
        % build the corresponding selection string
        str_sel = '';
        for i=1:(s-1)
            str_sel = [str_sel,int2str(ind(i))];
            str_sel = [str_sel,','];
        end;
        str_sel = [str_sel,':'];
        if( s~=d ) str_sel = [str_sel,',']; end;
        for i=(s+1):d
            str_sel = [str_sel,int2str(ind(i-1))];
            if( i~=d ) str_sel = [str_sel,',']; end;
        end;
        
        % put the vector in v
        str_stock = ['v=y(',str_sel, ');'];
        eval(str_stock);
        
        % reshape v into a 1D row vector
        v = reshape(v,1,length(v));
        % disp(str);
        
        % perform the transform
        str = ['v=', func_name, '(v'];
        if( isempty(func_arg) )
            str = [str, ');'];
        else
            str = [str, ',', func_arg, ');'];
        end
        % disp(str);
        eval(str);
        
        % reshape the result
        str = ['sz = size(y(', str_sel, '));'];
        eval(str);
        v = reshape(v,sz);
        
        % at last, assign the values        
        str = ['y(', str_sel, ')=v;'];
%        disp(str);
        eval(str);
    end
end