function y = perform_dct2_transform(x,dir,dim)

% perform_dct2_transform - perform the DCT2 transform
%
%   y = perform_dct2_transform(x);
%
%   Copyright (c) 2007

if nargin<2
    dir=1;
end
if nargin<3
    dim = 1;
    if size(x,1)>1 & size(x,2)>1
        dim = 2;
    end
end

n = size(x,1);
if dim==1
    % perform 1D DCT
    p = prod(size(x))/n;
    if 0
        %% implementation with size x 4 FFT
        y = zeros(4*n,p);
        y(2:2:2*n,:) = x(:,:); 
        y(4*n:-2:2*n+1,:) = x(:,:);
        y = fft(y,[],1);
        y = real( y(1:n,:) ) / sqrt(2*n);
        y(1,:) = y(1,:)/sqrt(2);
    else
        %% implementation with size x 2 FFT
        w = exp( -1i*pi/(2*n)*(0:n-1)');
        y = cat(1, x(:,:),zeros(n,p) );
        y = fft(y, [], 1); 
        y = y(1:n,:);
        y = sqrt(2/n) * real( y.*repmat( w, [1 p] ) ); 
        y(1,:) = y(1,:)/sqrt(2);
    end 
    y = reshape(y,size(x));
    return;
end

%% apply on 1st dimension
s = size(x); x = x(:,:); 
y = perform_dct2_transform(x,dir,1);
y = permute(y,[2 1 3:10]);
y = perform_dct2_transform(y,dir,1);
y = permute(y,[2 1 3:10]);

return;

if size(x,2)==1 || dim==1
    % 1D transform
    y = zeros(4*n,size(x,2),size(x,3),size(x,4));
    y(2:2:2*n,:,:,:) = x; y(4*n:-2:2*n+1,:,:,:) = x;
    y = fft(y,[],1);
    y = real( y(1:n,:,:,:) ) / sqrt(2*n);
    y(:,1,:,:) = y(:,1,:,:)/sqrt(2);
else
    if dir==1
        if  0
        y = x;
        y(1:end/2,:,:,:) = x(1:2:end,:,:,:);
        y(end:-1:end/2+1,:,:,:) = x(2:2:end,:,:,:);
        y = fft(y,[],1);
        y = real( y.*repmat( exp(-1i*pi/(2*n)*(0:n-1)'), [1,size(x,2),size(x,3),size(x,4)]  ) );
        end
        
        % transform on x
        y = zeros(4*n,size(x,2),size(x,3),size(x,4));
        y(2:2:2*n,:,:,:) = x; y(4*n:-2:2*n+1,:,:,:) = x;
        y = fft(y,[],1);
        y = real( y(1:n,:,:,:) );
        % transform on y
        z = zeros(size(x,1),4*n,size(x,3),size(x,4));
        z(:,2:2:2*n,:,:) = y; z(:,4*n:-2:2*n+1,:,:) = y;
        z = fft(z,[],2);
        y = real( z(:,1:n,:,:) ) / (2*n);
        % y(:,1,:,:) = y(:,1,:,:)/sqrt(2);
        % y(1,:,:,:) = y(1,:,:,:)/sqrt(2);
    else
        error('Not yet implemented (DCT3)');

    end
end