function [img1] = perform_bilateral_filtering(img,winsize,sigma)

% Bilateral Filtering(img,winsize,sigma)
% Input -> Image img
% -> winsize: spatial filter width
% -> sigma for intensity diff gaussain filter
% -> sigma for spatial filter = winsize/6
% Output -> Filtered Image
% Author: Amit Agrawal, 2004

% should be odd
winsize = round( (winsize-1)/2 )*2 + 1;

[H,W] = size(img);
%Gaussian spatial filter
g_filter = fspecial('gaussian',winsize,winsize/6);
padnum = (winsize-1)/2;
A = padarray(img,[padnum padnum],'replicate','both');
img1 = zeros(size(img));
for jj = padnum+1:(padnum+1+H-1)
    for kk = padnum+1:(padnum+1+W-1)
        % Get a local neighborhood
        imgwin = A(jj-padnum:jj+padnum,kk-padnum:kk+padnum);
        % Find weights according to intensity diffs
        Wwin = exp(-abs(imgwin - imgwin(padnum+1,padnum+1))/sigma^2);
        % Find composite filter
        newW = Wwin.*g_filter;
        t = sum(sum(newW));
        if(t>0)
            newW = newW/t;
        end
        img1(jj-padnum,kk-padnum) = sum(sum(imgwin.*newW));
    end
end