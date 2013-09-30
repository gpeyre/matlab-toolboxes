function y = callback_fft(x,dir,options)

% callback_fft - callback for sparsity with FFT
%
%   y = callback_fft(x,dir,options);
%
%   Works in 1D and 2D. Orthogonal transforms.
%
%   Copyright (c) 2008 Gabriel Peyre

options.null = 0;

%% Detect dimension
if size(x,1)==1 || size(x,2)==1
    ndims = 1;
else
    ndims = 2;
end
ndims = getoptions(options, 'ndims', ndims);
isreal = getoptions(options, 'isreal', +1);
remove_high_freq = getoptions(options, 'remove_high_freq', 0);



%% Transform
if ndims==1
    %% 1D
    n = length(x);
    if dir==1
        y = ifft(x) * sqrt(n);
        if isreal
            y = real(y);
        end
    else
        y = fft(x) / sqrt(n);
    end
else
    n = size(x,1); p = size(x,2);
    %% 2D
    if dir==1      
        if remove_high_freq>0
            x = remove_freq(x,remove_high_freq,dir);          
        end
        y = ifft2(x) * sqrt(n*p);
        if isreal
            y = real(y);
        end
    else
        y = fft2(x) / sqrt(n*p);
        if remove_high_freq>0
            y = remove_freq(y,remove_high_freq,dir);          
        end
    end
end

%%
function y = remove_freq(y,s,dir)


%y = fftshift(y);
% y(1:s,:) = 0; y(:,1:s) = 0;
% y(end-s+2:end,:) = 0; y(:,end-s+2:end) = 0;
% y = fftshift(y);

% remove only corners
y([end-s+2:end 1:s], end/2-s:end/2+s) = 0;
y(end/2-s:end/2+s, [end-s+2:end 1:s]) = 0;
y(end/2-s:end/2+s, end/2-s:end/2+s) = 0;


return;


global Qshuf; 
global Ishuf;
if size(Qshuf)~=3*(2*s-1)^2
    A = zeros(size(y));
    n = size(y,1);
    s1 = n/2-s+1:n/2+s-1;
    s2 = [n-s+2:n 1:s];
    A(s2, s1) = 1;
    A(s1, s2) = 1;
    A(s1,s1) = 1;
    Ishuf = find(A==1);
    rand('state', 123456);
    [Qshuf,R] = qr(rand(length(Ishuf)));
end
if dir==1
    y(Ishuf) = Qshuf*y(Ishuf);
else
    y(Ishuf) = Qshuf'*y(Ishuf);
end


return;
