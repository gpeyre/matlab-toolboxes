function y = perform_waveatoms_transform(x,dir, options)

% perform_waveatoms_transform - interface to WaveAtom transform
%
%   y = perform_waveatoms_transform(x,dir, options);
%
%   The waveatom toolbox can be downloaded from
%       http://www.waveatom.org/
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
issym = getoptions(options, 'issym', 0);
force_real = getoptions(options, 'force_real', 1);

if dir==1
    d = nb_dims(x);
else
    if iscell(x)
        d = nb_dims(x{1});
    else
        d = nb_dims(x);
    end
end

pat = 'p'; tp = 'ortho';
if d==1
    if dir==1
        if issym
            y = fatom1sym(x,pat,[1 1]);
        else
            y = fwa1(x,pat,tp);
        end
    else
        if issym
            y = iatom1sym(x,pat,[1 1]);
        else
            y = iwa1(x,pat,tp);
        end
    end
else
    if dir==1
        if issym
            y = fwa2sym(x,pat,[1 1]);
        else
            y = fwa2(x,pat,tp);
        end
    else
        if issym
            y = iwa2sym(x,pat,[1 1]);
        else
            y = iwa2(x,pat,tp);
        end
    end
end

if force_real && dir==-1
    y = real(y);
end


function c = fwa2(x,pat,tp)
% fwa2 - 2D forward wave atom transform
% -----------------
% INPUT
% --
% x is a real N-by-N matrix. N is a power of 2.
% --
% pat specifies the type of frequency partition which satsifies
% parabolic scaling relationship. pat can either be 'p' or 'q'.
% --
% tp is the type of tranform.
% 	'ortho': orthobasis
% 	'directional': real-valued frame with single oscillation direction
% 	'complex': complex-valued frame
% -----------------
% OUTPUT
% --
% c is a cell array which contains the wave atom coefficients. If
% tp=='ortho', then c{j}{m1,m2}(n1,n2) is the coefficient at scale j,
% frequency index (m1,m2) and spatial index (n1,n2). If
% tp=='directional', then c{j,d}{m1,m2}(n1,n2) with d=1,2 are the
% coefficients at scale j, frequency index (m1,m2) and spatial index
% (n1,n2). If tp=='complex', then c{j,d}{m1,m2)(n1,n2) with d=1,2,3,4
% are the coefficients at scale j, frequency index (m1,m2) and spatial
% index (n1,n2).
% -----------------
% Written by Lexing Ying and Laurent Demanet, 2007

if( ismember(tp, {'ortho','directional','complex'})==0 | ismember(pat, {'p','q','u'})==0 )    error('wrong');  end

if(strcmp(tp, 'ortho')==1)
    %---------------------------------------------------------
    N = size(x,1);
    H = N/2;
    lst = freq_pat(H,pat);
    %------------------
    f = fft2(x) / sqrt(prod(size(x)));
    A = N;
    c = cell(length(lst),1);
    %------------------
    for s=1:length(lst)
        nw = length(lst{s});
        c{s} = cell(nw,nw);
        for I=0:nw-1
            for J=0:nw-1
                if(lst{s}(I+1)==0 & lst{s}(J+1)==0)
                    c{s}{I+1,J+1} = [];
                else
                    B = 2^(s-1);
                    D = 2*B;
                    Ict = I*B;      Jct = J*B; %starting position in freq
                    if(mod(I,2)==0)
                        Ifm = Ict-2/3*B;        Ito = Ict+4/3*B;
                    else
                        Ifm = Ict-1/3*B;        Ito = Ict+5/3*B;
                    end
                    if(mod(J,2)==0)
                        Jfm = Jct-2/3*B;        Jto = Jct+4/3*B;
                    else
                        Jfm = Jct-1/3*B;        Jto = Jct+5/3*B;
                    end
                    res = zeros(D,D);
                    for id=0:1
                        if(id==0)
                            Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
                        else
                            Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
                        end
                        for jd=0:1
                            if(jd==0)
                                Jdx = [ceil(Jfm):floor(Jto)];      Jcf = kf_rt(Jdx/B*pi, J);
                            else
                                Jdx = [ceil(-Jto):floor(-Jfm)];      Jcf = kf_lf(Jdx/B*pi, J);
                            end
                            res(mod(Idx,D)+1,mod(Jdx,D)+1) = res(mod(Idx,D)+1,mod(Jdx,D)+1) + conj( Icf.'*Jcf ) .* f(mod(Idx,A)+1,mod(Jdx,A)+1);
                        end
                    end
                    c{s}{I+1,J+1} = ifft2(res) * sqrt(prod(size(res)));
                end
            end
        end
    end

elseif(strcmp(tp, 'directional')==1)
    %---------------------------------------------------------
    N = size(x,1);
    H = N/2;
    lst = freq_pat(H,pat);
    %------------------
    f = fft2(x) / sqrt(prod(size(x)));
    A = N;
    c1 = cell(length(lst),1);
    c2 = cell(length(lst),1);
    %------------------
    for s=1:length(lst)
        nw = length(lst{s});
        c1{s} = cell(nw,nw);
        c2{s} = cell(nw,nw);
        for I=0:nw-1
            for J=0:nw-1
                if(lst{s}(I+1)==0 & lst{s}(J+1)==0)
                    c1{s}{I+1,J+1} = [];
                    c2{s}{I+1,J+1} = [];
                else
                    B = 2^(s-1);
                    D = 2*B;
                    res = zeros(D,D);
                    Ict = I*B;      Jct = J*B; %starting position in freq
                    if(mod(I,2)==0)
                        Ifm = Ict-2/3*B;        Ito = Ict+4/3*B;
                    else
                        Ifm = Ict-1/3*B;        Ito = Ict+5/3*B;
                    end
                    if(mod(J,2)==0)
                        Jfm = Jct-2/3*B;        Jto = Jct+4/3*B;
                    else
                        Jfm = Jct-1/3*B;        Jto = Jct+5/3*B;
                    end

                    res = zeros(D,D);
                    Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
                    Jdx = [ceil(Jfm):floor(Jto)];      Jcf = kf_rt(Jdx/B*pi, J);
                    res(mod(Idx,D)+1,mod(Jdx,D)+1) = res(mod(Idx,D)+1,mod(Jdx,D)+1) + conj( Icf.'*Jcf ) .* f(mod(Idx,A)+1,mod(Jdx,A)+1);
                    Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
                    Jdx = [ceil(-Jto):floor(-Jfm)];      Jcf = kf_lf(Jdx/B*pi, J);
                    res(mod(Idx,D)+1,mod(Jdx,D)+1) = res(mod(Idx,D)+1,mod(Jdx,D)+1) + conj( Icf.'*Jcf ) .* f(mod(Idx,A)+1,mod(Jdx,A)+1);
                    c1{s}{I+1,J+1} = ifft2(res) * sqrt(prod(size(res)));

                    res = zeros(D,D);
                    Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
                    Jdx = [ceil(-Jto):floor(-Jfm)];      Jcf = kf_lf(Jdx/B*pi, J);
                    res(mod(Idx,D)+1,mod(Jdx,D)+1) = res(mod(Idx,D)+1,mod(Jdx,D)+1) + conj( Icf.'*Jcf ) .* f(mod(Idx,A)+1,mod(Jdx,A)+1);
                    Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
                    Jdx = [ceil(Jfm):floor(Jto)];      Jcf = kf_rt(Jdx/B*pi, J);
                    res(mod(Idx,D)+1,mod(Jdx,D)+1) = res(mod(Idx,D)+1,mod(Jdx,D)+1) + conj( Icf.'*Jcf ) .* f(mod(Idx,A)+1,mod(Jdx,A)+1);
                    c2{s}{I+1,J+1} = ifft2(res) * sqrt(prod(size(res)));
                end
            end
        end
    end
    c = [c1 c2];

elseif(strcmp(tp, 'complex')==1)
    %---------------------------------------------------------
    N = size(x,1);
    H = N/2;
    lst = freq_pat(H,pat);
    %------------------
    f = fft2(x) / sqrt(prod(size(x)));
    A = N;
    c1 = cell(length(lst),1);
    c2 = cell(length(lst),1);
    c3 = cell(length(lst),1);
    c4 = cell(length(lst),1);
    %------------------
    for s=1:length(lst)
        nw = length(lst{s});
        c1{s} = cell(nw,nw);
        c2{s} = cell(nw,nw);
        c3{s} = cell(nw,nw);
        c4{s} = cell(nw,nw);
        for I=0:nw-1
            for J=0:nw-1
                if(lst{s}(I+1)==0 & lst{s}(J+1)==0)
                    c1{s}{I+1,J+1} = [];
                    c2{s}{I+1,J+1} = [];
                    c3{s}{I+1,J+1} = [];
                    c4{s}{I+1,J+1} = [];
                else
                    B = 2^(s-1);
                    D = 2*B;
                    res = zeros(D,D);
                    Ict = I*B;      Jct = J*B; %starting position in freq
                    if(mod(I,2)==0)
                        Ifm = Ict-2/3*B;        Ito = Ict+4/3*B;
                    else
                        Ifm = Ict-1/3*B;        Ito = Ict+5/3*B;
                    end
                    if(mod(J,2)==0)
                        Jfm = Jct-2/3*B;        Jto = Jct+4/3*B;
                    else
                        Jfm = Jct-1/3*B;        Jto = Jct+5/3*B;
                    end

                    res = zeros(D,D);
                    Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
                    Jdx = [ceil(Jfm):floor(Jto)];      Jcf = kf_rt(Jdx/B*pi, J);
                    res(mod(Idx,D)+1,mod(Jdx,D)+1) = res(mod(Idx,D)+1,mod(Jdx,D)+1) + conj( Icf.'*Jcf ) .* f(mod(Idx,A)+1,mod(Jdx,A)+1);
                    c1{s}{I+1,J+1} = ifft2(res) * sqrt(prod(size(res)));

                    res = zeros(D,D);
                    Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
                    Jdx = [ceil(-Jto):floor(-Jfm)];      Jcf = kf_lf(Jdx/B*pi, J);
                    res(mod(Idx,D)+1,mod(Jdx,D)+1) = res(mod(Idx,D)+1,mod(Jdx,D)+1) + conj( Icf.'*Jcf ) .* f(mod(Idx,A)+1,mod(Jdx,A)+1);
                    c2{s}{I+1,J+1} = ifft2(res) * sqrt(prod(size(res)));

                    res = zeros(D,D);
                    Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
                    Jdx = [ceil(-Jto):floor(-Jfm)];      Jcf = kf_lf(Jdx/B*pi, J);
                    res(mod(Idx,D)+1,mod(Jdx,D)+1) = res(mod(Idx,D)+1,mod(Jdx,D)+1) + conj( Icf.'*Jcf ) .* f(mod(Idx,A)+1,mod(Jdx,A)+1);
                    c3{s}{I+1,J+1} = ifft2(res) * sqrt(prod(size(res)));

                    res = zeros(D,D);
                    Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
                    Jdx = [ceil(Jfm):floor(Jto)];      Jcf = kf_rt(Jdx/B*pi, J);
                    res(mod(Idx,D)+1,mod(Jdx,D)+1) = res(mod(Idx,D)+1,mod(Jdx,D)+1) + conj( Icf.'*Jcf ) .* f(mod(Idx,A)+1,mod(Jdx,A)+1);
                    c4{s}{I+1,J+1} = ifft2(res) * sqrt(prod(size(res)));
                end
            end
        end
    end
    c = [c1 c2 c3 c4];
end


function x = iwa2(c,pat,tp)
% iwa2 - 2D inverse wave atom transform
% -----------------
% INPUT
% --
% c is a cell array which contains the wave atom coefficients. If
% tp=='ortho', then c{j}{m1,m2}(n1,n2) is the coefficient at scale j,
% frequency index (m1,m2) and spatial index (n1,n2). If
% tp=='directional', then c{j,d}{m1,m2}(n1,n2) with d=1,2 are the
% coefficients at scale j, frequency index (m1,m2) and spatial index
% (n1,n2). If tp=='complex', then c{j,d}{m1,m2)(n1,n2) with d=1,2,3,4
% are the coefficients at scale j, frequency index (m1,m2) and spatial
% index (n1,n2).
% --
% pat specifies the type of frequency partition which satsifies
% parabolic scaling relationship. pat can either be 'p' or 'q'.
% --
% tp is the type of tranform.
% 	'ortho': orthobasis
% 	'directional': real-valued frame with single oscillation direction
% 	'complex': complex-valued frame
% -----------------
% OUTPUT
% --
% x is a real N-by-N matrix. N is a power of 2.
% -----------------
% Written by Lexing Ying and Laurent Demanet, 2007

if( ismember(tp, {'ortho','directional','complex'})==0 | ismember(pat, {'p','q','u'})==0 )    error('wrong');  end

if(strcmp(tp, 'ortho')==1)
    %---------------------------------------------------------
    T = 0;
    for s=1:length(c)
        nw = length(c{s});
        for I=1:nw
            for J=1:nw
                T = T + prod(size(c{s}{I,J}));
            end
        end
    end
    N = sqrt(T);
    H = N/2;
    lst = freq_pat(H,pat);
    A = N;
    f = zeros(A,A);
    %------------------
    for s=1:length(lst)
        nw = length(lst{s});
        for I=0:nw-1
            for J=0:nw-1
                if(~isempty(c{s}{I+1,J+1}))
                    B = 2^(s-1);
                    D = 2*B;
                    Ict = I*B;      Jct = J*B; %starting position in freq
                    if(mod(I,2)==0)
                        Ifm = Ict-2/3*B;        Ito = Ict+4/3*B;
                    else
                        Ifm = Ict-1/3*B;        Ito = Ict+5/3*B;
                    end
                    if(mod(J,2)==0)
                        Jfm = Jct-2/3*B;        Jto = Jct+4/3*B;
                    else
                        Jfm = Jct-1/3*B;        Jto = Jct+5/3*B;
                    end
                    res = fft2(c{s}{I+1,J+1}) / sqrt(prod(size(c{s}{I+1,J+1}))); %res = zeros(D,D);
                    for id=0:1
                        if(id==0)
                            Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
                        else
                            Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
                        end
                        for jd=0:1
                            if(jd==0)
                                Jdx = [ceil(Jfm):floor(Jto)];      Jcf = kf_rt(Jdx/B*pi, J);
                            else
                                Jdx = [ceil(-Jto):floor(-Jfm)];      Jcf = kf_lf(Jdx/B*pi, J);
                            end
                            f(mod(Idx,A)+1,mod(Jdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1) + ( Icf.'*Jcf ) .* res(mod(Idx,D)+1,mod(Jdx,D)+1);
                        end
                    end
                end
            end
        end
    end
    %------------------
    x = ifft2(f) * sqrt(prod(size(f)));

elseif(strcmp(tp, 'directional')==1)
    %---------------------------------------------------------
    c1 = c(:,1);
    c2 = c(:,2);

    T = 0;
    for s=1:length(c1)
        nw = length(c1{s});
        for I=1:nw
            for J=1:nw
                T = T + prod(size(c1{s}{I,J}));
            end
        end
    end
    N = sqrt(T);
    H = N/2;
    lst = freq_pat(H,pat);
    A = N;
    f = zeros(A,A);
    %------------------
    for s=1:length(lst)
        nw = length(lst{s});
        for I=0:nw-1
            for J=0:nw-1
                if(~isempty(c1{s}{I+1,J+1}))
                    B = 2^(s-1);
                    D = 2*B;
                    Ict = I*B;      Jct = J*B; %starting position in freq
                    if(mod(I,2)==0)
                        Ifm = Ict-2/3*B;        Ito = Ict+4/3*B;
                    else
                        Ifm = Ict-1/3*B;        Ito = Ict+5/3*B;
                    end
                    if(mod(J,2)==0)
                        Jfm = Jct-2/3*B;        Jto = Jct+4/3*B;
                    else
                        Jfm = Jct-1/3*B;        Jto = Jct+5/3*B;
                    end

                    res = fft2(c1{s}{I+1,J+1}) / sqrt(prod(size(c1{s}{I+1,J+1}))); %res = zeros(D,D);
                    Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
                    Jdx = [ceil(Jfm):floor(Jto)];      Jcf = kf_rt(Jdx/B*pi, J);
                    f(mod(Idx,A)+1,mod(Jdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1) + ( Icf.'*Jcf ) .* res(mod(Idx,D)+1,mod(Jdx,D)+1);
                    Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
                    Jdx = [ceil(-Jto):floor(-Jfm)];      Jcf = kf_lf(Jdx/B*pi, J);
                    f(mod(Idx,A)+1,mod(Jdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1) + ( Icf.'*Jcf ) .* res(mod(Idx,D)+1,mod(Jdx,D)+1);

                    res = fft2(c2{s}{I+1,J+1}) / sqrt(prod(size(c2{s}{I+1,J+1}))); %res = zeros(D,D);
                    Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
                    Jdx = [ceil(-Jto):floor(-Jfm)];      Jcf = kf_lf(Jdx/B*pi, J);
                    f(mod(Idx,A)+1,mod(Jdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1) + ( Icf.'*Jcf ) .* res(mod(Idx,D)+1,mod(Jdx,D)+1);
                    Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
                    Jdx = [ceil(Jfm):floor(Jto)];      Jcf = kf_rt(Jdx/B*pi, J);
                    f(mod(Idx,A)+1,mod(Jdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1) + ( Icf.'*Jcf ) .* res(mod(Idx,D)+1,mod(Jdx,D)+1);
                end
            end
        end
    end
    x = ifft2(f) * sqrt(prod(size(f)));

elseif(strcmp(tp, 'complex')==1)
    %---------------------------------------------------------
    c1 = c(:,1);
    c2 = c(:,2);
    c3 = c(:,3);
    c4 = c(:,4);

    T = 0;
    for s=1:length(c1)
        nw = length(c1{s});
        for I=1:nw
            for J=1:nw
                T = T + prod(size(c1{s}{I,J}));
            end
        end
    end
    N = sqrt(T);
    H = N/2;
    lst = freq_pat(H,pat);
    A = N;
    f = zeros(A,A);
    %------------------
    for s=1:length(lst)
        nw = length(lst{s});
        for I=0:nw-1
            for J=0:nw-1
                if(~isempty(c1{s}{I+1,J+1}))
                    B = 2^(s-1);
                    D = 2*B;
                    Ict = I*B;      Jct = J*B; %starting position in freq
                    if(mod(I,2)==0)
                        Ifm = Ict-2/3*B;        Ito = Ict+4/3*B;
                    else
                        Ifm = Ict-1/3*B;        Ito = Ict+5/3*B;
                    end
                    if(mod(J,2)==0)
                        Jfm = Jct-2/3*B;        Jto = Jct+4/3*B;
                    else
                        Jfm = Jct-1/3*B;        Jto = Jct+5/3*B;
                    end

                    res = fft2(c1{s}{I+1,J+1}) / sqrt(prod(size(c1{s}{I+1,J+1}))); %res = zeros(D,D);
                    Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
                    Jdx = [ceil(Jfm):floor(Jto)];      Jcf = kf_rt(Jdx/B*pi, J);
                    f(mod(Idx,A)+1,mod(Jdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1) + ( Icf.'*Jcf ) .* res(mod(Idx,D)+1,mod(Jdx,D)+1);

                    res = fft2(c2{s}{I+1,J+1}) / sqrt(prod(size(c2{s}{I+1,J+1}))); %res = zeros(D,D);
                    Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
                    Jdx = [ceil(-Jto):floor(-Jfm)];      Jcf = kf_lf(Jdx/B*pi, J);
                    f(mod(Idx,A)+1,mod(Jdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1) + ( Icf.'*Jcf ) .* res(mod(Idx,D)+1,mod(Jdx,D)+1);

                    res = fft2(c3{s}{I+1,J+1}) / sqrt(prod(size(c3{s}{I+1,J+1}))); %res = zeros(D,D);
                    Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);
                    Jdx = [ceil(-Jto):floor(-Jfm)];      Jcf = kf_lf(Jdx/B*pi, J);
                    f(mod(Idx,A)+1,mod(Jdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1) + ( Icf.'*Jcf ) .* res(mod(Idx,D)+1,mod(Jdx,D)+1);

                    res = fft2(c4{s}{I+1,J+1}) / sqrt(prod(size(c4{s}{I+1,J+1}))); %res = zeros(D,D);
                    Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);
                    Jdx = [ceil(Jfm):floor(Jto)];      Jcf = kf_rt(Jdx/B*pi, J);
                    f(mod(Idx,A)+1,mod(Jdx,A)+1) = f(mod(Idx,A)+1,mod(Jdx,A)+1) + ( Icf.'*Jcf ) .* res(mod(Idx,D)+1,mod(Jdx,D)+1);
                end
            end
        end
    end
    x = ifft2(f) * sqrt(prod(size(f)));
end


function lst = freq_pat(H,pat)
% freq_pat.m - generate frequency partition
%
% input:
%   H:    half length of the frequency span, H needs to be greater than 16
%   pat:  type of frequency partition which satsifies parabolic scaling relationship
%         equal to 'p' or 'q'
%
% output:
%   lst:  data representing the partition of frequency
%
% Written by Lexing Ying and Laurent Demanet, 2006

if(H<16)    error('H needs to be at least 16');  end

if(    pat=='p') %parabolic scaling by lexing
    len = ceil(log2(H)/2)+1;
    lst = cell(len,1);
    lst{1} = [1 1];    lst{2} = [0 1];    lst{3} = [0 1 1 1];
    cnt = 16;
    idx = 3;    rad = 4;
    while(cnt<H)
        old = cnt;
        lst{idx} = [lst{idx} 1 1];
        cnt = cnt + 2*rad;
        idx = idx+1;
        rad = 2*rad;
        trg = min(4*old,H);
        lst{idx} = [zeros(1,cnt/rad), ones(1,(trg-cnt)/rad)];
        cnt = trg;
    end
elseif(pat=='q') %parabolic scaling by laurent
    len = floor(log2(H)/2)+1;
    lst = cell(len,1);
    lst{1} = [1 1];    lst{2} = [0 1 1 1];
    cnt = 8;
    idx = 2;      rad = 2;
    while(cnt<H)
        old = cnt;
        lst{idx} = [lst{idx} 1 1];
        cnt = cnt + 2*rad;
        idx = idx+1;
        rad = 2*rad;
        trg = min(4*old,H);
        lst{idx} = [zeros(1,cnt/rad), ones(1,(trg-cnt)/rad)];
        cnt = trg;
    end
elseif(pat=='u')
    len = ceil(log2(H)/2)+1; %uniform partitioning
    lst = cell(len,1);
    B = 2^(len-1);
    lst{end} = ones(1,H/B);
else
    error('wrong pat');
end

function r = kf_rt(w,n)
% kf_rt.m - right bump
%
% Written by Lexing Ying and Laurent Demanet, 2006

r = zeros(size(w));
an = pi/2*(n+1/2);
en = (-1)^n;
en1 = (-1)^(n+1);

r = exp(-i*w/2) .* ( exp(i*an)*g(en*(w-pi*(n+1/2))) );

function r = g(w)
% g.m - 'g' function Villemoes's construction
%
% Written by Lexing Ying and Laurent Demanet, 2006
r = zeros(size(w));
gd = w<5*pi/6 & w>-7*pi/6;
r(gd) = abs(sf(w(gd)-3*pi/2));

%----------------------------------------------------------------------
function r = sf(w)
r = zeros(size(w));

aw = abs(w);
r(aw<=2*pi/3) = 0;

gd = aw>=2*pi/3 & aw<=4*pi/3;
r(gd) = 1/sqrt(2) * hf(w(gd)/2+pi);

gd = aw>=4*pi/3 & aw<=8*pi/3;
r(gd) = 1/sqrt(2) * hf(w(gd)/4);

r(aw>8*pi/2) = 0;

function r = hf(w)
w = mod(w+pi,2*pi) - pi;
r = zeros(size(w));
w = abs(w);

r = sqrt(2) * cos(pi/2 * beta(3*w/pi-1));
r(w<=pi/3) = sqrt(2);
r(w>=2*pi/3) = 0;

function r = beta(x)
r = x.^4 .*(35-84*x + 70*x.^2 -20*x.^3);

function r = kf_lf(w,n)
% kf_lf.m - left bump
%
% Written by Lexing Ying and Laurent Demanet, 2006

r = zeros(size(w));
an = pi/2*(n+1/2);
en = (-1)^n;
en1 = (-1)^(n+1);

r = exp(-i*w/2) .* ( exp(-i*an)*g(en1*(w+pi*(n+1/2))) );

function res = fwa2sym(x,pat,tp)
% fwa2sym - 2D forward wave atom transform (symmetric version)
% -----------------
% INPUT
% --
% x is a real N-by-N matrix. N is a power of 2.
% --
% pat specifies the type of frequency partition which satsifies
% parabolic scaling relationship. pat can either be 'p' or 'q'.
% --
% tp is the type of tranform.
% 	'ortho': orthobasis
% 	'directional': real-valued frame with single oscillation direction
% 	'complex': complex-valued frame
% -----------------
% OUTPUT
% --
% res is an array containing all the wave atom coefficients. If
% tp=='ortho', then res is of size N-by-N. If tp=='
% If tp=='directional', then res is of size N-by-N-by-2. If tp=='complex',
% then res is of size N-by-N-by-4.
% -----------------
% Written by Lexing Ying and Laurent Demanet, 2007

if( ismember(tp, {'ortho','directional','complex'})==0 | ismember(pat, {'p','q','u'})==0 )    error('wrong');  end

call = fwa2(x,pat,tp);
N = size(x,1);
res = zeros(N,N,size(call,2));

for i=1:size(call,2)
    c = call(:,i);
    y = zeros(N,N);
    for s=1:length(c)
        D = 2^s;
        nw = length(c{s});
        for I=0:nw-1
            for J=0:nw-1
                if(~isempty(c{s}{I+1,J+1}))
                    y( I*D+[1:D], J*D+[1:D] ) = c{s}{I+1,J+1};
                end
            end
        end
    end
    res(:,:,i) = y;
end




function x = iwa2sym(res,pat,tp)
% iwa2sym - 2D inverse wave atom transform (symmetric version)
% -----------------
% INPUT
% --
% res is an array containing all the wave atom coefficients. If
% tp=='ortho', then res is of size N-by-N. If tp=='
% If tp=='directional', then res is of size N-by-N-by-2. If tp=='complex',
% then res is of size N-by-N-by-4.
% --
% pat specifies the type of frequency partition which satsifies
% parabolic scaling relationship. pat can either be 'p' or 'q'.
% --
% tp is the type of tranform.
% 	'ortho': orthobasis
% 	'directional': real-valued frame with single oscillation direction
% 	'complex': complex-valued frame
% -----------------
% OUTPUT
% --
% x is a real N-by-N matrix. N is a power of 2.
% -----------------
% Written by Lexing Ying and Laurent Demanet, 2007


if( ismember(tp, {'ortho','directional','complex'})==0 | ismember(pat, {'p','q','u'})==0 )    error('wrong');  end

N = size(res,1);
H = N/2;
lst = freq_pat(H,pat);
red = size(res,3);

call = cell(length(lst),red);
for i=1:red
    y = res(:,:,i);
    c = cell(length(lst),1);
    for s=1:length(lst)
        B = 2^(s-1);    D = 2*B;
        nw = length(lst{s});
        c{s} = cell(nw,nw);
        for I=0:nw-1
            for J=0:nw-1
                if(lst{s}(I+1)==0 & lst{s}(J+1)==0)
                    c{s}{I+1,J+1} = [];
                else
                    c{s}{I+1,J+1} = y( I*D+[1:D], J*D+[1:D] );
                end
            end
        end
    end
    call(:,i) = c;
end
x = iwa2(call,pat,tp);



