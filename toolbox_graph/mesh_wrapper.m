function dist = mesh_wrapper(vertex1,face1,vertex2,face2);

% mesh_wrapper - calls the MESH software and parse the result
%
%   dist = metro_wrapper(vertex1,face1,vertex2,face2);
%
%   Copyright (c) 2004 Gabriel Peyré

% write the two off files
write_smf( 'tmp1.smf', vertex1,face1 );
write_smf( 'tmp2.smf', vertex2,face2 );

nb_samples = size(face1,1);

% launch MESH tool
str = '!mesh -q -s -t tmp1.smf tmp2.smf > tmp.mesh';
tic;
eval(str);
disp( sprintf('Mesh took %.1fs.', toc) );


% parsing of the file
fid = fopen('tmp.mesh','rt');
if( fid==-1 )
    error('Can''t open the file.');
    return;
end

min = [];
min_bb = [];
max = [];
max_bb = [];
mean = [];
mean_bb = [];
rms = [];
rms_bb = [];

str = 1;
while( 1 )
    str = fgets(fid);   % -1 if eof
    if str==-1
        dist.min = min;
        dist.min_bb = min_bb;
        dist.max = max;
        dist.max_bb = max_bb;
        dist.mean = mean;
        dist.mean_bb = mean_bb;
        dist.rms = rms;
        dist.rms_bb = rms_bb;
        return;
    end
    [a,str] = strtok(str);
    if strcmp(a, 'Min:')
        [a,str] = strtok(str); min = [min, str2double(a)];
        [a,str] = strtok(str); min_bb = [min_bb, str2double(a)];
    end
    if strcmp(a, 'Max:')
        [a,str] = strtok(str); max = [max, str2double(a)];
        [a,str] = strtok(str); max_bb = [max_bb, str2double(a)];
    end
    if strcmp(a, 'Mean:')
        [a,str] = strtok(str); mean = [mean, str2double(a)];
        [a,str] = strtok(str); mean_bb = [mean_bb, str2double(a)];
    end
    if strcmp(a, 'RMS:')
        [a,str] = strtok(str); rms = [rms, str2double(a)];
        [a,str] = strtok(str); rms_bb = [rms_bb, str2double(a)];
    end
end


fclose(fid);

% delete the results
!del tmp*.smf
!del tmp.mesh