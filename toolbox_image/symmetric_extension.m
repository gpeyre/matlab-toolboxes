function M_padded = symmetric_extension(M,k)

% symmetric_extension - perform a symmetric extension of the signal.
%
%   M_padded = symmetric_extension(M,k);
%
% M can be 1D or 2D array
% If size(M)=[n,p], the result is of size [n+2*k,p+2*k]
%
%   Copyright (c) 2004 Gabriel PeyrŽ

if size(M,3)>1
    M_padded = zeros( size(M,1)+2*k, size(M,2)+2*k, size(M,3) );
    for i=1:size(M,3)
        M_padded(:,:,i) = symmetric_extension(M(:,:,i),k);
    end
    return;
end

n1 = size(M,1);
n2 = size(M,2);

if nb_dims(M)==1
    M = M(:);
    M_padded = [ M(k:-1:1); M; M(end:-1:end-k+1) ];
elseif nb_dims(M)==2
    M_padded = zeros(n1+2*k,n2+2*k);
    M_padded(k+1:end-k,k+1:end-k) = M;
    % extension
    M_padded(1:k,:) = M_padded(2*k:-1:k+1,:);
    M_padded(end-k+1:end,:) = M_padded(end-k:-1:end-2*k+1,:);
    M_padded(:,1:k) = M_padded(:,2*k:-1:k+1);
    M_padded(:,end-k+1:end) = M_padded(:,end-k:-1:end-2*k+1);
else
    error('Only supported for array of dimension less than 2.')
end


function k = nb_dims(x)

if size(x,1)==1 || size(x,2)==1
    k = 1;
else
    k=2;
end