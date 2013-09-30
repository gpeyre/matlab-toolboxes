function [y,nbr_bits] = perform_spiht_coding(x,options)

% perform_spiht_coding - SPIHT coding of wavelet coefficients
%
% Coding : 
%   options.Jmin = ??;      % minimum scale of the transform
%   options.nb_bits = ??;   % target number of bits
%   [stream,nbr_bits] = perform_spiht_coding(MW,options);
% Decoding : 
%   MW = perform_spiht_coding(stream);
%
% This is a simple wrapper of the code of Jing Tian
% <scuteejtian at hotmail.com>

options.null = 0;

if isfield(options, 'arithmetic_coding')
    arithmetic_coding = options.arithmetic_coding;
else
    arithmetic_coding = 1;
end

if size(x,1)>1 && size(x,2)>1
    if isfield(options, 'Jmin')
        Jmin = options.Jmin;
    else
        warning('You should provide options.Jmin');
        Jmin = 4;
    end
    if isfield(options, 'nb_bits')
        nb_bits = options.nb_bits;
    else
        nb_bits = floor( prod(size(x))*0.1 );
    end
    %-----------   Coding   ----------------
    Jmax = log2(size(x,1))-1;
    level = Jmax-Jmin+1;
    y = func_SPIHT_Enc(x, nb_bits, level); y = y(:);
    if arithmetic_coding
        % remove trailing 2
        I = find(y==2); y(I) = [];
        % the 3 first entry are [size nb_bitsplanes level]
        [z,nbr_bits] = perform_arithmetic_coding(y(4:end), +1, options);
        y = [y(1:3); z];
        nbr_bits = nbr_bits + 16; % approx 16 bits
    else
        nbr_bits = nb_bits;
    end
else
    %-----------   Decoding   ----------------
    if arithmetic_coding
        x = x(:);
        % the 3 first entry are [size nb_bitsplanes level]
        z = perform_arithmetic_coding(x(4:end), -1, options);
        x = [x(1:3); z];
        nbr_bits = -1;
    end
    y = func_SPIHT_Dec(x(:)');
end

function out = func_SPIHT_Enc(m, max_bits, level)
% Matlab implementation of SPIHT (without Arithmatic coding stage)
%
% Encoder
%
% input:    m : input image in wavelet domain
%           max_bits : maximum bits can be used
%           level : wavelet decomposition level
%
% output:   out : bit stream
%
% Jing Tian
% Contact me : scuteejtian@hotmail.com
% This program is part of my undergraduate project in GuangZhou, P. R. China.
% April - July 1999


%-----------   Initialization  -----------------
bitctr = 0;
out = 2*ones(1,max_bits - 14);
n_max = floor(log2(abs(max(max(m)'))));
Bits_Header = 0;
Bits_LSP = 0;
Bits_LIP = 0;
Bits_LIS = 0;

%-----------   output bit stream header   ----------------
% image size, number of bit plane, wavelet decomposition level should be
% written as bit stream header.
out(1,[1 2 3]) = [size(m,1) n_max level]; bitctr = bitctr + 24;
index = 4;
Bits_Header = Bits_Header + 24;

%-----------   Initialize LIP, LSP, LIS   ----------------
temp = [];
bandsize = 2.^(log2(size(m, 1)) - level + 1);
temp1 = 1 : bandsize;
for i = 1 : bandsize
    temp = [temp; temp1];
end
LIP(:, 1) = temp(:);
temp = temp';
LIP(:, 2) = temp(:);
LIS(:, 1) = LIP(:, 1);
LIS(:, 2) = LIP(:, 2);
LIS(:, 3) = zeros(length(LIP(:, 1)), 1);
pstart = 1;
pend = bandsize / 2;
for i = 1 : bandsize / 2
    LIS(pstart : pend, :) = [];
    pdel = pend - pstart + 1;
    pstart = pstart + bandsize - pdel;
    pend = pend + bandsize - pdel;
end
LSP = [];

n = n_max;

%-----------   coding   ----------------
while(bitctr < max_bits)
        
    % Sorting Pass
    LIPtemp = LIP; temp = 0;
    for i = 1:size(LIPtemp,1)
        temp = temp+1;
        if (bitctr + 1) >= max_bits
            if (bitctr < max_bits)
                out(length(out))=[];
            end
            return
        end
        if abs(m(LIPtemp(i,1),LIPtemp(i,2))) >= 2^n % 1: positive; 0: negative
            out(index) = 1; bitctr = bitctr + 1;
            index = index +1; Bits_LIP = Bits_LIP + 1;
            sgn = m(LIPtemp(i,1),LIPtemp(i,2))>=0;
            out(index) = sgn; bitctr = bitctr + 1;
            index = index +1; Bits_LIP = Bits_LIP + 1;
            LSP = [LSP; LIPtemp(i,:)];
            LIP(temp,:) = []; temp = temp - 1;
        else
            out(index) = 0; bitctr = bitctr + 1;
            index = index +1;
            Bits_LIP = Bits_LIP + 1;
        end
    end
    
    LIStemp = LIS; temp = 0; i = 1;
    while ( i <= size(LIStemp,1))
        temp = temp + 1;
        if LIStemp(i,3) == 0
            if bitctr >= max_bits
                return
            end
            max_d = func_MyDescendant(LIStemp(i,1),LIStemp(i,2),LIStemp(i,3),m);
            if max_d >= 2^n
                out(index) = 1; bitctr = bitctr + 1;
                index = index +1; Bits_LIS = Bits_LIS + 1;
                x = LIStemp(i,1); y = LIStemp(i,2);
                
                if (bitctr + 1) >= max_bits
                    if (bitctr < max_bits)
                        out(length(out))=[];
                    end
                    return
                end
                if abs(m(2*x-1,2*y-1)) >= 2^n
                    LSP = [LSP; 2*x-1 2*y-1];
                    out(index) = 1; bitctr = bitctr + 1;
                    index = index +1; Bits_LIS = Bits_LIS + 1;
                    sgn = m(2*x-1,2*y-1)>=0;
                    out(index) = sgn; bitctr = bitctr + 1;
                    index = index +1; Bits_LIS = Bits_LIS + 1;
                else
                    out(index) = 0; bitctr = bitctr + 1;
                    index = index +1; Bits_LIS = Bits_LIS + 1;
                    LIP = [LIP; 2*x-1 2*y-1];
                end
                
                if (bitctr + 1) >= max_bits
                    if (bitctr < max_bits)
                        out(length(out))=[];
                    end
                    return
                end
                if abs(m(2*x-1,2*y)) >= 2^n
                    LSP = [LSP; 2*x-1 2*y];
                    out(index) = 1; bitctr = bitctr + 1;
                    index = index +1; Bits_LIS = Bits_LIS + 1;
                    sgn = m(2*x-1,2*y)>=0;
                    out(index) = sgn; bitctr = bitctr + 1;
                    index = index +1; Bits_LIS = Bits_LIS + 1;
                else
                    out(index) = 0; bitctr = bitctr + 1;
                    index = index +1; Bits_LIS = Bits_LIS + 1;
                    LIP = [LIP; 2*x-1 2*y];
                end
                
                if (bitctr + 1) >= max_bits
                    if (bitctr < max_bits)
                        out(length(out))=[];
                    end
                    return
                end
                if abs(m(2*x,2*y-1)) >= 2^n
                    LSP = [LSP; 2*x 2*y-1];
                    out(index) = 1; bitctr = bitctr + 1;
                    index = index +1; Bits_LIS = Bits_LIS + 1;
                    sgn = m(2*x,2*y-1)>=0;
                    out(index) = sgn; bitctr = bitctr + 1;
                    index = index +1; Bits_LIS = Bits_LIS + 1;
                else
                    out(index) = 0; bitctr = bitctr + 1;
                    index = index +1; Bits_LIS = Bits_LIS + 1;
                    LIP = [LIP; 2*x 2*y-1];
                end
                
                if (bitctr + 1) >= max_bits
                    if (bitctr < max_bits)
                        out(length(out))=[];
                    end
                    return
                end
                if abs(m(2*x,2*y)) >= 2^n
                    LSP = [LSP; 2*x 2*y];
                    out(index) = 1; bitctr = bitctr + 1;
                    index = index +1; Bits_LIS = Bits_LIS + 1;
                    sgn = m(2*x,2*y)>=0;
                    out(index) = sgn; bitctr = bitctr + 1;
                    index = index +1; Bits_LIS = Bits_LIS + 1;
                else
                    out(index) = 0; bitctr = bitctr + 1;
                    index = index +1; Bits_LIS = Bits_LIS + 1;
                    LIP = [LIP; 2*x 2*y];
                end
                
                if ((2*(2*x)-1) < size(m) & (2*(2*y)-1) < size(m))
                    LIS = [LIS; LIStemp(i,1) LIStemp(i,2) 1];
                    LIStemp = [LIStemp; LIStemp(i,1) LIStemp(i,2) 1];
                end
                LIS(temp,:) = []; temp = temp-1;
                
            else
                out(index) = 0; bitctr = bitctr + 1;
                index = index +1; Bits_LIS = Bits_LIS + 1;
            end
        else
            if bitctr >= max_bits
                return
            end
            max_d = func_MyDescendant(LIStemp(i,1),LIStemp(i,2),LIStemp(i,3),m);
            if max_d >= 2^n
                out(index) = 1; bitctr = bitctr + 1;
                index = index +1;
                x = LIStemp(i,1); y = LIStemp(i,2);
                LIS = [LIS; 2*x-1 2*y-1 0; 2*x-1 2*y 0; 2*x 2*y-1 0; 2*x 2*y 0];
                LIStemp = [LIStemp; 2*x-1 2*y-1 0; 2*x-1 2*y 0; 2*x 2*y-1 0; 2*x 2*y 0];
                LIS(temp,:) = []; temp = temp - 1;
            else
                out(index) = 0; bitctr = bitctr + 1;
                index = index +1; Bits_LIS = Bits_LIS + 1;
            end
        end
        i = i+1;
    end
    
    % Refinement Pass
    temp = 1;
    value = floor(abs(2^(n_max-n+1)*m(LSP(temp,1),LSP(temp,2))));
    while (value >= 2^(n_max+2) & (temp <= size(LSP,1)))
        if bitctr >= max_bits
            return
        end
        s = bitget(value,n_max+2);
        out(index) = s; bitctr = bitctr + 1;
        index = index +1; Bits_LSP = Bits_LSP + 1;
        temp = temp + 1;
        if temp <= size(LSP,1)
            value = floor(abs(2^(n_max-n+1)*m(LSP(temp,1),LSP(temp,2))));
        end
    end
    
    n = n - 1;
end




function m = func_SPIHT_Dec(in)
% Matlab implementation of SPIHT (without Arithmatic coding stage)
%
% Decoder
%
% input:    in : bit stream
%
% output:   m : reconstructed image in wavelet domain
%
% Jing Tian
% Contact me : scuteejtian@hotmail.com
% This program is part of my undergraduate project in GuangZhou, P. R. China.
% April - July 1999

%-----------   Initialization  -----------------
% image size, number of bit plane, wavelet decomposition level should be
% written as bit stream header.
m = zeros(in(1,1));
n_max = in(1,2);
level = in(1,3);
ctr = 4;
 
%-----------   Initialize LIP, LSP, LIS   ----------------
temp = [];
bandsize = 2.^(log2(in(1,1)) - level + 1);
temp1 = 1 : bandsize;
for i = 1 : bandsize
    temp = [temp; temp1];
end
LIP(:, 1) = temp(:);
temp = temp';
LIP(:, 2) = temp(:);

LIS(:, 1) = LIP(:, 1);
LIS(:, 2) = LIP(:, 2);
LIS(:, 3) = zeros(length(LIP(:, 1)), 1);
pstart = 1;
pend = bandsize / 2;
for i = 1 : bandsize / 2
    LIS(pstart : pend, :) = [];
    pdel = pend - pstart + 1;
    pstart = pstart + bandsize - pdel;
    pend = pend + bandsize - pdel;
end
LSP = [];

%-----------   coding   ----------------
n = n_max;
while (ctr <= size(in,2))
    
    %Sorting Pass
    LIPtemp = LIP; temp = 0;
    for i = 1:size(LIPtemp,1)
        temp = temp+1;
        if ctr > size(in,2)
            return
        end
        if in(1,ctr) == 1
            ctr = ctr + 1;
            if in(1,ctr) > 0
                m(LIPtemp(i,1),LIPtemp(i,2)) = 2^n + 2^(n-1);  
            else
                m(LIPtemp(i,1),LIPtemp(i,2)) = -2^n  - 2^(n-1); 
            end
            LSP = [LSP; LIPtemp(i,:)];
            LIP(temp,:) = []; temp = temp - 1;
        end
        ctr = ctr + 1;
    end
    
    LIStemp = LIS; temp = 0; i = 1;
    while ( i <= size(LIStemp,1))
        temp = temp + 1;
        if ctr > size(in,2)
            return
        end
        if LIStemp(i,3) == 0
            if in(1,ctr) == 1 
                ctr = ctr + 1;
                x = LIStemp(i,1); y = LIStemp(i,2);
                
                if ctr > size(in,2)
                    return
                end
                if in(1,ctr) == 1
                    LSP = [LSP; 2*x-1 2*y-1];
                    ctr = ctr + 1;
                    if in(1,ctr) == 1
                        m(2*x-1,2*y-1) = 2^n + 2^(n-1); 
                    else
                        m(2*x-1,2*y-1) = -2^n  - 2^(n-1); 
                    end
                    ctr = ctr + 1;
                else
                    LIP = [LIP; 2*x-1 2*y-1];
                    ctr = ctr + 1;
                end
                
                if ctr > size(in,2)
                    return
                end
                if in(1,ctr) == 1
                    ctr = ctr + 1;
                    LSP = [LSP; 2*x-1 2*y];
                    if in(1,ctr) == 1;
                        m(2*x-1,2*y) = 2^n + 2^(n-1); 
                    else
                        m(2*x-1,2*y) = -2^n  - 2^(n-1); 
                    end
                    ctr = ctr + 1;
                else
                    LIP = [LIP; 2*x-1 2*y];
                    ctr = ctr + 1;
                end
                
                if ctr > size(in,2)
                    return
                end
                if in(1,ctr) == 1
                    ctr = ctr + 1;
                    LSP = [LSP; 2*x 2*y-1];
                    if in(1,ctr) == 1
                        m(2*x,2*y-1) = 2^n + 2^(n-1); 
                    else
                        m(2*x,2*y-1) = -2^n  - 2^(n-1);
                    end
                    ctr = ctr + 1;
                else
                    LIP = [LIP; 2*x 2*y-1];
                    ctr = ctr + 1;
                end
                
                if ctr > size(in,2)
                    return
                end
                if in(1,ctr) == 1
                    ctr = ctr + 1;
                    LSP = [LSP; 2*x 2*y];
                    if in(1,ctr) == 1
                        m(2*x,2*y) = 2^n + 2^(n-1); 
                    else
                        m(2*x,2*y) = -2^n  - 2^(n-1); 
                    end
                    ctr = ctr + 1;
                else
                    LIP = [LIP; 2*x 2*y];
                    ctr = ctr + 1;
                end
                
                if ((2*(2*x)-1) < size(m) & (2*(2*y)-1) < size(m))
                    LIS = [LIS; LIStemp(i,1) LIStemp(i,2) 1];
                    LIStemp = [LIStemp; LIStemp(i,1) LIStemp(i,2) 1];
                end
                LIS(temp,:) = []; temp = temp-1;
                
            else
                ctr = ctr + 1;
            end
        else
            if in(1,ctr) == 1
                x = LIStemp(i,1); y = LIStemp(i,2);
                LIS = [LIS; 2*x-1 2*y-1 0; 2*x-1 2*y 0; 2*x 2*y-1 0; 2*x 2*y 0];
                LIStemp = [LIStemp; 2*x-1 2*y-1 0; 2*x-1 2*y 0; 2*x 2*y-1 0; 2*x 2*y 0];
                LIS(temp,:) = []; temp = temp - 1;
            end
            ctr = ctr + 1;
        end
        i = i+1;
    end
    
    % Refinement Pass
    temp = 1;
    value = m(LSP(temp,1), LSP(temp,2));
    while (abs(value) >= 2^(n+1) & (temp <= size(LSP,1)))
        if ctr > size(in,2)
            return
        end

        value = value + ((-1)^(in(1,ctr) + 1)) * (2^(n-1))*sign(m(LSP(temp,1),LSP(temp,2))); 
        m(LSP(temp,1),LSP(temp,2)) = value;
        ctr = ctr + 1;
        temp = temp + 1;    
        if temp <= size(LSP,1)
            value = m(LSP(temp,1),LSP(temp,2));
        end
    end
    
    n = n-1;
end




function value = func_MyDescendant(i, j, type, m)
% Matlab implementation of SPIHT (without Arithmatic coding stage)
%
% Find the descendant with largest absolute value of pixel (i,j)
%
% input:    i : row coordinate
%           j : column coordinate
%           type : type of descendant
%           m : whole image
%
% output:   value : largest absolute value
%
% Jing Tian
% Contact me : scuteejtian@hotmail.com
% This program is part of my undergraduate project in GuangZhou, P. R. China.
% April - July 1999

s = size(m,1);

S = [];

index = 0; a = 0; b = 0;

while ((2*i-1)<s & (2*j-1)<s)
    a = i-1; b = j-1;

    mind = [2*(a+1)-1:2*(a+2^index)];
    nind = [2*(b+1)-1:2*(b+2^index)];
    
    
    chk = mind <= s;
    len = sum(chk);
    if len < length(mind)
        mind(len+1:length(mind)) = [];
    end
    
    
    chk = nind <= s;
    len = sum(chk);
    if len < length(nind)
        nind(len+1:length(nind)) = [];
    end
    
    S = [S reshape(m(mind,nind),1,[])];
    
    index = index + 1;
    i = 2*a+1; j = 2*b+1;
end

if type == 1
    S(:,1:4) = [];; 
end

value = max(abs(S));