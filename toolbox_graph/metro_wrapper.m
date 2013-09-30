function dist = metro_wrapper(vertex1,face1, vertex2,face2)

% metro_wrapper - calls the metro software and parse the result
%
%   dist = metro_wrapper(vertex1,face1,vertex2,face2);
%
%   Copyright (c) 2004 Gabriel Peyré

% write the two off files
write_off( 'tmp1.off', vertex1,face1 );
write_off( 'tmp2.off', vertex2,face2 );

write_smf( 'tmp1.smf', vertex1,face1 );
write_smf( 'tmp2.smf', vertex2,face2 );

nb_samples = size(face1,1);

% launch metro tool
str = ['!metro tmp1.off tmp2.off -n' num2str(nb_samples) ' > tmp.metro'];
tic;
eval(str);
disp( sprintf('Metro took %.1fs.', toc) );


% parsing of the file
fid = fopen('tmp.metro','rt');
if( fid==-1 )
    error('Can''t open the file.');
    return;
end

str = 1;
i = 1;
while( str ~= -1)
    str = fgets(fid);   % -1 if eof
    if( ~isempty(findstr(str,'distances:')) )
        str = fgets(fid);
        [a,str] = strtok(str);
        [a,str] = strtok(str);
        [a,str] = strtok(str); max(i) = str2double(a);
        [a,str] = strtok(str); max_bb(i) = str2double(a(2:end));
        str = fgets(fid);
        [a,str] = strtok(str);
        [a,str] = strtok(str);
        [a,str] = strtok(str); mean(i) = str2double(a);
        str = fgets(fid);
        [a,str] = strtok(str);
        [a,str] = strtok(str);
        [a,str] = strtok(str); rms(i) = str2double(a);
        i = i+1;
    end
end

dist.max = max;
dist.max_bb = max_bb;
dist.mean = mean;
dist.rms = rms;

fclose(fid);

% delete the results
!del tmp*.off
!del tmp.metro