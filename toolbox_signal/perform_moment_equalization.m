function x = perform_moment_equalization(x,y,numdim, options)

% perform_kurtosis_equalization - equalize moments of order 1,2,3,4.
%
%   x = perform_moment_equalization(x,y,numdim,options);
%
%   (numdim=1 by default).
%
%   Equalizes the mean, variance, skewness and kurtosis.
%   Set options.xx=0 to avoid equalizing one of the moment, where
%       xx='skewness' or 'kurtosis'.
%
%   For 2D arrays, operates along columns unless you specify
%       x = perform_kurtosis_equalization(x,y,2);
%
% Copyright (c) 2006 Gabriel Peyr?

%   mem=mean2(chm);
%   sk2 = mean2((chm-mem).^3)/mean2((chm-mem).^2).^(3/2);

options.null = 0;
if nargin<3
    numdim = 1;
end

if size(x,1)>1 && size(x,2)>1
    if numdim>1
        x = x'; y = y';
    end
    p = size(x,2);
    for i=1:p
        if p>100
            progressbar(i,p);
        end
        x(:,i) = perform_moment_equalization( x(:,i),y(:,i),1,options );
    end
    if numdim>1
        x = x';
    end
    return;
end

dokurt = 1;
doskew = 1;
if isfield(options, 'kurtosis')
    dokurt = options.kurtosis;
end
if isfield(options, 'skewness')
    doskew = options.skewness;
end

niter = 1;


if std(x(:))<eps && std(y(:))>eps
    % regenerate random noise ...
    x = randn(size(x)) * std(y(:)) + mean(y(:));
end

if std(y(:))>eps
    sk = skew2(y);
    k = kurt2(y);
    for i=1:niter
        if dokurt
            [x, snrk] = modkurt(x,k);
        end
        if doskew
            [x, snrk] = modskew(x,sk);
        end
    end
end
% correct mean and variance
if std(x(:))>1e-9
    x = (x-mean(x(:))) * std(y(:))/std(x(:)) + mean(y(:));
end

function [chm, snrk] = modkurt(ch,k,p);

% Modify the kurtosis in one step, by moving in gradient direction until
% reaching the desired kurtosis value. 
% It does not affect the mean nor the variance, but it affects the skewness.
% This operation is not an orthogonal projection, but the projection angle is
% near pi/2 when k is close to the original kurtosis, which is a realistic assumption
% when doing iterative projections in a pyramid, for example (small corrections
% to the channels' statistics). 
%
% [chm, snrk] = modkurt(ch,k,p);
%	ch: channel
%	k: desired kurtosis (k=M4/M2^2)	
%       p [OPTIONAL]:   mixing proportion between k0 and k
%                       it imposes (1-p)*k0 + p*k,
%                       being k0 the current kurtosis.
%                       DEFAULT: p = 1;


% Javier Portilla,  Oct.12/97, NYU

Warn = 0;  % Set to 1 if you want to see warning messages
if ~exist('p'),
  p = 1;
end

me=mean2(ch);
ch=ch-me;

% Compute the moments

m=zeros(12,1);
for n=2:12,
	m(n)=mean2(ch.^n);
end

% The original kurtosis

k0=m(4)/m(2)^2;
snrk = snr(k, k-k0);
if snrk > 60,
	chm = ch+me;
	return
end
k = k0*(1-p) + k*p;

% Some auxiliar variables

a=m(4)/m(2);

% Coeficients of the numerator (A*lam^4+B*lam^3+C*lam^2+D*lam+E)

A=m(12)-4*a*m(10)-4*m(3)*m(9)+6*a^2*m(8)+12*a*m(3)*m(7)+6*m(3)^2*m(6)-...
	4*a^3*m(6)-12*a^2*m(3)*m(5)+a^4*m(4)-12*a*m(3)^2*m(4)+...
	4*a^3*m(3)^2+6*a^2*m(3)^2*m(2)-3*m(3)^4;
B=4*(m(10)-3*a*m(8)-3*m(3)*m(7)+3*a^2*m(6)+6*a*m(3)*m(5)+3*m(3)^2*m(4)-...
	a^3*m(4)-3*a^2*m(3)^2-3*m(4)*m(3)^2);
C=6*(m(8)-2*a*m(6)-2*m(3)*m(5)+a^2*m(4)+2*a*m(3)^2+m(3)^2*m(2));
D=4*(m(6)-a^2*m(2)-m(3)^2);
E=m(4);

% Define the coefficients of the denominator (F*lam^2+G)^2

F=D/4;
G=m(2);

% test
test = 0;

if test,

        grd = ch.^3 - a*ch - m(3);
        lam = -0.001:0.00001:0.001;
        k = (A*lam.^4+B*lam.^3+C*lam.^2+D*lam+E)./...
                (F*lam.^2 + G).^2;
        for lam = -0.001:0.00001:0.001,
                n = lam*100000+101;
                chp = ch + lam*grd;
                k2(n) = mean2(chp.^4)/mean2(chp.^2)^2;
                %k2(n) = mean2(chp.^4);
        end
        lam = -0.001:0.00001:0.001;
        snr(k2, k-k2)

end % test

% Now I compute its derivative with respect to lambda
% (only the roots of derivative = 0 )

d(1) = B*F;
d(2) = 2*C*F - 4*A*G;
d(3) = 4*F*D -3*B*G - D*F;
d(4) = 4*F*E - 2*C*G;
d(5) = -D*G;

mMlambda = roots(d);

tg = imag(mMlambda)./real(mMlambda);
mMlambda = mMlambda(find(abs(tg)<1e-6));
lNeg = mMlambda(find(mMlambda<0));
if length(lNeg)==0,
        lNeg = -1/eps;
end
lPos = mMlambda(find(mMlambda>=0));
if length(lPos)==0,
        lPos = 1/eps;
end
lmi = max(lNeg);
lma = min(lPos);

lam = [lmi lma];
mMnewKt = polyval([A B C D E],lam)./(polyval([F 0 G],lam)).^2;
kmin = min(mMnewKt);
kmax = max(mMnewKt);

% Given a desired kurtosis, solves for lambda

if k<=kmin
        lam = lmi;
        if Warn
        warning('Saturating (down) kurtosis!');
        kmin
        end
elseif k>=kmax
        lam = lma;
        if Warn
        warning('Saturating (up) kurtosis!');
        kmax
        end
else

% Coeficients of the algebraic equation

c0 = E - k*G^2;
c1 = D;
c2 = C - 2*k*F*G;
c3 = B;
c4 = A - k*F^2;

% Solves the equation

r=roots([c4 c3 c2 c1 c0]);

% Chose the real solution with minimum absolute value with the rigth sign

denom = real(r);
denom(abs(denom)<eps) = 1;
tg = imag(r)./denom;

%lambda = real(r(find(abs(tg)<1e-6)));
lambda = real(r(find(abs(tg)==0)));
if length(lambda)>0,
	lam = lambda(find(abs(lambda)==min(abs(lambda))));
	lam = lam(1);
else
	lam = 0;
end

end % if ... else


% Modify the channel

chm=ch+lam*(ch.^3-a*ch-m(3));	% adjust the kurtosis
chm=chm*sqrt(m(2)/mean2(chm.^2));	% adjust the variance
chm=chm+me;				% adjust the mean

% Check the result
%k2=mean2((chm-me).^4)/(mean2((chm-me).^2))^2;
%SNR=snr(k,k-k2)








 


function [chm, snrk] = modskew(ch,sk,p);

% Adjust the sample skewness of a vector/matrix, using gradient projection,
% without affecting its sample mean and variance.
%
% This operation is not an orthogonal projection, but the projection angle is
% near pi/2 when sk is close to the original skewness, which is a realistic
% assumption when doing iterative projections in a pyramid, for example
% (small corrections to the channels' statistics).
%
% 	[xm, snrk] = modskew(x,sk,p);
%		sk: new skweness
%       	p [OPTIONAL]:   mixing proportion between sk0 and sk 
%               	        it imposes (1-p)*sk0 + p*sk,
%                      		being sk0 the current skewness.
%                       	DEFAULT: p = 1;

%
% JPM. 2/98, IODV, CSIC
% 4/00, CNS, NYU

Warn = 0;  % Set to 1 if you want to see warning messages
if ~exist('p'), 
  p = 1;
end

N=prod(size(ch));	% number of samples
me=mean2(ch);
ch=ch-me;

for n=2:6,
	m(n)=mean2(ch.^n);
end

sd=sqrt(m(2));	% standard deviation
s=m(3)/sd^3;	% original skewness
snrk = snr(sk, sk-s); 
sk = s*(1-p) + sk*p;

% Define the coefficients of the numerator (A*lam^3+B*lam^2+C*lam+D)

A=m(6)-3*sd*s*m(5)+3*sd^2*(s^2-1)*m(4)+sd^6*(2+3*s^2-s^4);
B=3*(m(5)-2*sd*s*m(4)+sd^5*s^3);
C=3*(m(4)-sd^4*(1+s^2));
D=s*sd^3;

a(7)=A^2;
a(6)=2*A*B;
a(5)=B^2+2*A*C;
a(4)=2*(A*D+B*C);
a(3)=C^2+2*B*D;
a(2)=2*C*D;
a(1)=D^2;

% Define the coefficients of the denominator (A2+B2*lam^2)

A2=sd^2;
B2=m(4)-(1+s^2)*sd^4;

b=zeros(1,7);
b(7)=B2^3;
b(5)=3*A2*B2^2;
b(3)=3*A2^2*B2;
b(1)=A2^3;


if 0, % test

	lam = -2:0.02:2;
	S = (A*lam.^3+B*lam.^2+C*lam+D)./...
		sqrt(b(7)*lam.^6 + b(5)*lam.^4 + b(3)*lam.^2 + b(1));
%	grd = ch.^2 - m(2) - sd * s * ch;
%	for lam = -1:0.01:1,
%		n = lam*100+101;
%		chp = ch + lam*grd;
%		S2(n) = mean2(chp.^3)/abs(mean2(chp.^2))^(1.5);
%	end
	lam = -2:0.02:2;
figure(1);plot(lam,S);grid;drawnow
%	snr(S2, S-S2)

end % test

% Now I compute its derivative with respect to lambda

d(8) = B*b(7);
d(7) = 2*C*b(7) - A*b(5);
d(6) = 3*D*b(7);
d(5) = C*b(5) - 2*A*b(3);
d(4) = 2*D*b(5) - B*b(3);
d(3) = -3*A*b(1);
d(2) = D*b(3) - 2*B*b(1);
d(1) = -C*b(1);

d = d(8:-1:1);
mMlambda = roots(d);

denom = real(mMlambda);
[tmp,I] = find( min(abs(denom))<eps );
denom(I) = 1;
tg = imag(mMlambda)./denom;
mMlambda = real(mMlambda(find(abs(tg)<1e-6)));
lNeg = mMlambda(find(mMlambda<0));
if length(lNeg)==0,
	lNeg = -1/eps;
end
lPos = mMlambda(find(mMlambda>=0));
if length(lPos)==0,
        lPos = 1/eps;
end
lmi = max(lNeg);
lma = min(lPos);

lam = [lmi lma];
denom = (polyval(b(7:-1:1),lam)).^0.5;
if 0
[tmp,I] = find( min(abs(denom))<eps );
denom(I) = 1;
end
denom(abs(denom)<eps) = 1;

mMnewSt = polyval([A B C D],lam)./denom;
    
skmin = min(mMnewSt);
skmax = max(mMnewSt);


% Given a desired skewness, solves for lambda

if sk<=skmin
        lam = lmi;
        if Warn
            warning('Saturating (down) skewness!');
            skmin
        end
elseif sk>=skmax & Warn,
        lam = lma;
        if Warn
            warning('Saturating (up) skewness!');
            skmax
        end
else


% The equation is sum(c.*lam.^(0:6))=0

c=a-b*sk^2;

c=c(7:-1:1);

r=roots(c);

% Chose the real solution with minimum absolute value with the rigth sign
lam=-Inf;
co=0;
for n=1:min(6,length(r)),
    denom = real(r(n));
    [tmp,I] = find( min(abs(denom))<eps );
    denom(I) = 1;
	tg = imag(r(n))/denom;
	if (abs(tg)<1e-6)&(sign(real(r(n)))==sign(sk-s)),
		co=co+1;
		lam(co)=real(r(n));
	end
end
if min(abs(lam))==Inf
    if Warn
        display('Warning: Skew adjustment skipped!');
    end
	lam=0;
end

p=[A B C D];

if length(lam)>1,
	foo=sign(polyval(p,lam));
	if any(foo==0),
		lam = lam(find(foo==0));
	else
		lam = lam(find(foo==sign(sk)));		% rejects the symmetric solution
	end
	if length(lam)>0,
		lam=lam(find(abs(lam)==min(abs(lam))));	% the smallest that fix the skew
		lam=lam(1);
	else
		lam = 0;
	end
end
end % if else

% Modify the channel
chm=ch+lam*(ch.^2-sd^2-sd*s*ch);	% adjust the skewness
chm=chm*sqrt(m(2)/mean2(chm.^2));		% adjust the variance
chm=chm+me;				% adjust the mean
					% (These don't affect the skewness)
% Check the result
%mem=mean2(chm);
%sk2=mean2((chm-mem).^3)/mean2((chm-mem).^2).^(3/2);
%sk - sk2
%SNR=snr(sk,sk-sk2)





% S = SKEW2(MTX,MEAN,VAR)
%
% Sample skew (third moment divided by variance^3/2) of a matrix.
%  MEAN (optional) and VAR (optional) make the computation faster.

function res = skew2(mtx, mn, v)

if (exist('mn') ~= 1)
  mn =  mean2(mtx);
end

if (exist('v') ~= 1)
  v =  var2(mtx,mn);
end

if (isreal(mtx))
  res = mean(mean((mtx-mn).^3)) / (v^(3/2));
else
  res = mean(mean(real(mtx-mn).^3)) / (real(v)^(3/2)) + ...
      i * mean(mean(imag(mtx-mn).^3)) / (imag(v)^(3/2));
end


% K = KURT2(MTX,MEAN,VAR)
%
% Sample kurtosis (fourth moment divided by squared variance) 
% of a matrix.  Kurtosis of a Gaussian distribution is 3.
%  MEAN (optional) and VAR (optional) make the computation faster.

% Eero Simoncelli, 6/96.

function res = kurt2(mtx, mn, v)

if (exist('mn') ~= 1)
	mn =  mean(mean(mtx));
end

if (exist('v') ~= 1)
	v =  var2(mtx,mn);
end

if (isreal(mtx))
  res = mean(mean(abs(mtx-mn).^4)) / (v^2);
else
  res = mean(mean(real(mtx-mn).^4)) / (real(v)^2) + ...
      i*mean(mean(imag(mtx-mn).^4)) / (imag(v)^2);
end


% M = MEAN2(MTX)
%
% Sample mean of a matrix.

function res = mean2(mtx)

res = mean(mean(mtx));


% V = VAR2(MTX,MEAN)
%
% Sample variance of a matrix.
%  Passing MEAN (optional) makes the calculation faster.

function res = var2(mtx, mn)

if (exist('mn') ~= 1)
  mn =  mean2(mtx);
end

if (isreal(mtx))
  res = sum(sum(abs(mtx-mn).^2)) / max((prod(size(mtx)) - 1),1);
else
  res = sum(sum(real(mtx-mn).^2)) + i*sum(sum(imag(mtx-mn).^2));
  res = res  / max((prod(size(mtx)) - 1),1);
end


function X=snr(s,n);

% Compute the signal-to-noise ratio in dB
%  	X=SNR(signal,noise);
% (it does not subtract the means).

es=sum(sum(abs(s).^2));
en=sum(sum(abs(n).^2));
if en<eps
    X = en;
else
    X=10*log10(es/en);
end

