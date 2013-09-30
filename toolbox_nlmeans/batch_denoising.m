% batch for denoising

namelist = {'lena','barb','boat', 'peppers', 'mandrill', 'polygons_blurred'};
for iname=1:length(namelist);
    name = namelist{iname};
    disp(['----> Denoising ' name '.']);
    test_denoising;
end