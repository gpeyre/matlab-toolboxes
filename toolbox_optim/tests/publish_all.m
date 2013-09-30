%%
% Compile all HTML files.

a = dir('test_*.m');
opts.outputDir = 'html';

for i=1:length(a) 
    name = a(i).name(1:end-2);
	publish(name,opts);
end