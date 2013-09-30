% set the paths for all the tooboxes
%
%   Copyright (c) 2008 Gabriel Peyre

% return;

% base directory for the toolboxes
rep = 'Desktop/dev/matlab/matlabcompanion/toolbox/';
rep = '/Volumes/gpeyre 1/matlab/toolbox/';

    
% return;

if not(rep(1)=='/');

% retrieve the user name
s = pwd();
f = strfind(s, '/');
s = s(1:f(3));
rep = [s rep];

end

rep_old = pwd;
cd(rep);

A = dir('toolbox_*');
for i=1:length(A)
    name = A(i).name;
    if not(strcmp(name(end-2:end), 'zip'))
        path(path, [rep name]);
    end
end

% set default font size larger
set_thicklines(2);

cd(rep_old);