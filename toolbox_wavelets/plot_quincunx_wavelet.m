function MW = plot_quincunx_wavelet(MW, Jmin, options)

% plot_wavelet - plot 2D wavelet transform stored using Mallat's ordering.
%
%   plot_wavelet(MW, Jmin, options);
%
%   'MW' is the wavelet transform.
%   'Jmin' is the minimum scale of the transform.
%
%   Copyright (c) 2003 Gabriel Peyré

n = size(MW,1);
Jmax = log2(n)-1;
J = (Jmax-Jmin+1)*2;

bnd = 1; % boundary value

options.null = 0;
gamma = getoptions(options, 'gamma', 1);

% select some sub-images
options.transform = 'quincunx';
m = 2*(Jmax-Jmin+1)+1;
nb = ceil(m/2);
k = 0;
% clf;
for j=Jmax:-1:Jmin
    qq = 1:2;
    if  j==Jmin
        qq(end+1) = 0;
    end
    for q=qq
        k = k+1;
        [selx,sely] = compute_quadrant_selection(j,q, options);
        MWj = MW(selx,sely);
        if q==0
            MWj = rescale(MWj);
        else
            MWj = MWj/max(abs(MWj(:)));
            % gamma boosting of contrast
            MWj = abs(MWj.^gamma).*sign(MWj);
            MWj = (MWj+1)/2;
        end
        
        MWj(1,:) = bnd; MWj(end,:) = bnd; MWj(:,1) = bnd; MWj(:,end) = bnd;
        MW(selx,sely) = MWj;
    end
end
colormap gray(256);

imageplot(MW);



return;

Q2 = showdecomposition(MW,J);
warning off;
Q2=uint8(Q2/max(abs(Q2(:)))*128+128);
warning on;
imshow(Q2,[]);
axis image; axis off;


%
% Show decomposition with each subband rescaled
%

function Q2 = showdecomposition(Q,J);

mmax=max(abs(Q(:)));

x2=size(Q,1);
y2=size(Q,2);

x1=1;
y1=y2/2+1;

Q2=zeros(size(Q));

for iter=1:J,
    tmp=Q(x1:x2,y1:y2);
    tmp=tmp/max(abs(tmp(:)))*mmax;
    Q2(x1:x2,y1:y2)=tmp;
    if mod(iter,2)==1,
        Q2(x1:x2,y1)=mmax;
        y2=y1-1;
        y1=1;
        x1=x2/2+1;
    else
        Q2(x1,y1:y2)=mmax;
        x2=x1-1;
        x1=1;
        y1=y2/2+1;
    end;
end;

y1=1;

if mod(iter,2)==1,
    x1=1;
end;

tmp=Q(x1:x2,y1:y2);
tmp=min(Q2(:))+tmp/max(abs(tmp(:)))*(max(Q2(:))-min(Q2(:)));
Q2(x1:x2,y1:y2)=tmp;

