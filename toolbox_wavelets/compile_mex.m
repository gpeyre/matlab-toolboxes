% compile the mex files for the Wavelet toolbox.
%
%   Copyright (c) 2006 Gabriel Peyr??



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% JPEG2000
disp('---> Compiling Jpeg2000 mex files.');

rep = 'jp2k/src/';
strbase = 'mex ';
if ispc
    % windows special definitions
    strbase = [strbase ' -Dcompil_vcc -DWIN32 '];
end

% Compile jp2_class
if 0
files = {'jp2_codec.c' 'bio.c' 'dwt.c' 'j2k.c' 'mct.c' 'pi.c' 't2.c' 'tgt.c' ...
            'cio.c' 'fix.c' 'int.c' 'mqc.c' 't1.c' ...
            'tcd.c' 'image_jp2.c' 'liw_jp2_dll.c' 'liw_error.c'};
str = [strbase '-output jp2_class '];
for i=1:length(files)
    str = [str rep files{i} ' '];
end
eval(str);
end
    
% Compile jp2_class
files = {'jp2_codec.c' 'bio.c' 'dwt.c' 'j2k.c' 'mct.c' 'pi.c' 't2.c' ...
            'tgt.c' 'cio.c' 'fix.c' 'int.c' 'mqc.c' 't1.c' 'tcd.c' ...
            'image_jp2.c' 'liw_jp2_dll.c' 'liw_error.c'};
str = [strbase ' -DENCODE_ONLY -output perform_jp2k_encoding '];
for i=1:length(files) 
    str = [str rep files{i} ' '];
end
eval(str);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LIW atrou
disp('---> Compiling a trou wavelet transforms.');

mex cwpt2/cwpt2_btree.c cwpt2/btree.c  
mex cwpt2/cwpt2.c cwpt2/dyadic2.c cwpt2/atrous.c cwpt2/btree.c cwpt2/quad.c
mex cwpt2/cwpt2i.c cwpt2/dyadic2.c cwpt2/atrous.c cwpt2/btree.c cwpt2/quad.c
mex cwpt2/c3wpt.c cwpt2/dyadic2.c cwpt2/atrous.c cwpt2/btree.c cwpt2/quad.c
mex cwpt2/c3wpti.c cwpt2/dyadic2.c cwpt2/atrous.c cwpt2/btree.c cwpt2/quad.c

if 0
mex cwt.c dyadic2.c atrous.c btree.c quad.c
mex cwti.c dyadic2.c atrous.c btree.c quad.c
mex cwpt2_get_leaf_order.c btree.c
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Toolbox files
disp('---> Compiling wavelet transforms mex files.');

mex mex/perform_79_transform.cpp
mex mex/perform_haar_transform.cpp
mex mex/perform_lifting_transform.cpp


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Toolbox files
disp('---> Compiling matlabSpyrTool needed mex files.');
mex  mex/simoncelli/upConv.c mex/simoncelli/wrap.c mex/simoncelli/convolve.c mex/simoncelli/edges.c



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RWT
disp('---> Compiling rice wavelet toolbox mex files.');
mex rwt/mdwt.c rwt/mdwt_r.c
mex rwt/midwt.c rwt/midwt_r.c
mex rwt/mrdwt.c rwt/mrdwt_r.c
mex rwt/mirdwt.c rwt/mirdwt_r.c