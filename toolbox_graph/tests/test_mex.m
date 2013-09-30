% test for mex file

test_read_vrml = 1;

if test_read_vrml 
   filename = 'pawn.wrl';
   [vertex,face] = read_vrml(filename);
end