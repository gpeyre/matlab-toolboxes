function z = poly_root(P)

% poly_root - find the roots of a polynomial
%
%   z = root_low_poly(P);
%
%   Special code for degrees < 4.
%   This method can take a list of polynomials, 
%   i.e. P(:,i) is the coefficient of the ith polynomial
%   (fast parallel resolution).
%
%   Copyright (c) 2004 Gabriel Peyré

if size(P,1)==1 && size(P,2)~=1
    P = P';
end

if size(P,1)==2
    
    z = -P(1,:)./P(2,:);
    
elseif size(P,1)==3
    
    delta = P(2,:).^2 - 4*P(1,:).*P(3,:);
    z = [ (-P(2,:)-sqrt(delta))./( 2*P(3,:) ); (-P(2,:)+sqrt(delta))./( 2*P(3,:) )];
    
elseif size(P,1)==4
    
    % normalize polynomial
    for i=1:3
        P(i,:) = P(i,:)./P(4,:);
    end
    P(4,:) = 1;
    
    p = -P(3,:).*P(3,:)/9+P(2,:)/3;
    q = 2*P(3,:).^3/27-P(3,:).*P(2,:)/3+P(1,:);
    rr = poly_root([-p.^3; q; P(4,:)]);   % remember P(4,:) = 1
    u = ( abs(rr(1,:)).^(1/3) ) .* sign(rr(1,:)); 
    v = ( abs(rr(2,:)).^(1/3) ) .* sign(rr(2,:)); 
    z = [u+v-P(3,:)/3;
        -0.5*(u+v)-P(3,:)/3 + 0.5i*sqrt(3).*(u-v);
        -0.5*(u+v)-P(3,:)/3 - 0.5i*sqrt(3).*(u-v)];
    
elseif size(P,1)==5
    
    % to be done ...
    p = -3*P(4,:).*P(4,:)/8 + P(3,:);
    q = P(4,:).^3/8 - 0.5.*P(4,:).*P(3,:)+P(2,:);
    s = -3*P(4,:).^4/256 + P(4,:).*P(4,:).*P(3,:)/16 - 0.25*P(4,:).*P(2,:) + P(1,:);
    
    I1 = find( abs(q)>=eps );
    I2 = find( abs(q)>=eps );
    
    q1 = q(I1);
    p1 = p(I2);
    rr = poly_root( [-q1*q1/64; (p1*p1-4*s1)/16 0.5*p1; P(5,I1)] );
    k = find( abs(rr)>eps ); 
    A =2*sqrt(rr(k(1)));
    z1 = [ roots1( [0.5*(p+A*A+q/A); -1; P(5,:)] );
    poly_root([0.5*(p+A*A-q/A); A; P(5,:)])]-0.25*P(4,:);

    p2 = p(I2);
    s2 = s(I2);
    u =roots1([1,p,s]); 
    z2 =[sqrt(u); -sqrt(u)]-0.25*P(4);


    z = zeros(size(P,1)-1, size(P,2));
    z(:,I1) = z1;
    z(:,I1) = z2;

else
    
    z = zeros(size(P,1)-1, size(P,2));
    for i=1:size(P,2)
        z(:,i) = roots(P(:,i));
    end
    
end


% I = find( imag(z)~=0 );
% z(I) = NaN;