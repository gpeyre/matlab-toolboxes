name_list = {'dunes' 'chocolate' 'pointilliste'  'corral' 'mures' 'pasta'  'frenchfries' 'tomatoes' };

%  'olives' 

for iname = 1:length(name_list)
    name = name_list{iname};
    disp(['************ synthesizing ' name ' **************']);
    do_quilting = 1;
    k_schedule = 3;
    test_nl_synthesis;
    do_quilting = 0;
    k_schedule = [7 6 5 4 3 2];
    test_nl_synthesis;
end

return;


name_list = {'group-people'};
Tlist = [1e-9 0.15 0.06];

% 'olives', 

for iname = 1:length(name_list)
    name = name_list{iname};
    for do_patchwise = 1:-1:0
        for T = Tlist
            test_nl_synthesis;
        end
    end
end