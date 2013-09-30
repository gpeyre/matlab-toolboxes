function save_image(M, name, options)

% save_image - save an image or a set of images
%
%   save_image(M, name, options);
%
%   M and name can be cell arrays.
%
%   set options.clamp and options.rescale to clamp/rescale the image prior
%   to saving.
%
%   set options.format='bmp', 'jpg' ... to set the file format.
%
%   set options.base_str for the string that is appened to each filename.
%
%   Copyright (c) 2008 Gabriel Peyre


options.null = 0;
if iscell(M)
    if not(iscell(M)) || length(name)~=length(M)
        error('name should be a cell array as M');
    end
    for i=1:length(M)
        save_image(M{i}, name{i}, options);
    end
    return;
end

% trying to guess format
s = strfind(name, '.');
fmt = 'png';
if not(isempty(s))
    fmt = name(s(1)+1:end);
end
fmt = getoptions(options, 'format', fmt);

if iscell(fmt)
    for i=1:length(fmt)
        options.format = fmt{i};
        save_image(M, name, options);
    end
    return;
end

c = getoptions(options, 'clamp', []);
r = getoptions(options, 'rescale', []);
if c==1
    c = [0 1];
end
if r==1
    r = [0 1];
end
if not(isempty(c))
    M = clamp(M, c(1), c(2));
end
if not(isempty(r))
    M = rescale(M, r(1), r(2));
end
base_str = getoptions(options, 'base_str', '');

if strcmp(fmt, 'eps') 
    warning('eps is unsuported for now');
    return;
    % detection of image magick
    global image_magick_path;
    if isempty(image_magick_path)
        image_magick_path = '/Applications/ImageMagick-6.4.3/bin/';
    end
    if not(exist([image_magick_path 'convert']))
        warning('Image magick not installed, you need to specify global image_magick_path');
    end
    
    warning off;
    imwrite(M, [base_str name '.png'],'png');
    warning on;
    system([ image_magick_path 'convert ./' base_str name '.png ./' base_str name '.eps' ]);
    delete( [base_str name '.png'] );
    return;
end

warning off;
imwrite(M, [base_str name '.' fmt], fmt);
warning on;