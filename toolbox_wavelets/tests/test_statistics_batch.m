name_list = {'lena','barb','reptilskin','hair'};


for iname = 1:length(name_list)
    name = name_list{iname};
    wavtype = 'redun';
    test_statistics;
    wavtype = 'ortho';
    test_statistics;
end