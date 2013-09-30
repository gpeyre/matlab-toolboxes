function imden=den2d_block(x,qmf,thd,L,sigma)

if (exist('x')~=1), 	error('Provide input image'); 		end;
if (exist('qmf')~=1), 	qmf = MakeONFilter('Symmlet',6); 	end;
if (exist('thd')~=1), 	thd = 4.5;				end;

[nx ny] = size(x);
if nx~=ny
	disp('Matrix must be square');
	return;
end

n 	= nx;
J	= nextpow2(n);
if (exist('L')~=1), 	L   = floor(sqrt(2*J)); 		end; % Size of block.
scale 	= floor(log2(L));   % Coarsest scale L#2^Jc matches the size of the block.

% wc=FWT2_PO(x,scale,qmf);
wc = x;

for j = 0:J-1-scale
 HH = wc(2^(J-j-1)+1 : 2^(J-j), 2^(J-j-1)+1 : 2^(J-j));
 HH = block_partition(HH,L);
 %HHmask = sqrt(sum(HH.^2)) >= thd*L*sigma;
 HHmask = max(1-thd*log(n)*sigma^2./sum(HH.^2),0);
 HH = HH .* repmat(HHmask,L*L,1);
 HH = invblock_partition(HH,2^(J-j-1),L);
 wc(2^(J-j-1)+1 : 2^(J-j), 2^(J-j-1)+1 : 2^(J-j)) = HH;
 
 HL = wc(2^(J-j-1)+1 : 2^(J-j), 1 : 2^(J-j-1));
 HL = block_partition(HL,L);
 %HLmask = sqrt(sum(HL.^2)) >= thd*L*sigma;
 HLmask = max(1-thd*log(n)*sigma^2./sum(HL.^2),0);
 HL = HL .* repmat(HLmask,L*L,1);
 HL = invblock_partition(HL,2^(J-j-1),L);
 wc(2^(J-j-1)+1 : 2^(J-j), 1 : 2^(J-j-1)) = HL;
 
 LH = wc(1 : 2^(J-j-1), 2^(J-j-1)+1 : 2^(J-j));
 LH = block_partition(LH,L);
 %LHmask = sqrt(sum(LH.^2)) >= thd*L*sigma;
 LHmask = max(1-thd*log(n)*sigma^2./sum(LH.^2),0);
 LH = LH .* repmat(LHmask,L*L,1);
 LH = invblock_partition(LH,2^(J-j-1),L);
 wc(1 : 2^(J-j-1), 2^(J-j-1)+1 : 2^(J-j)) = LH;
end

imden = IWT2_PO(wc,scale,qmf);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = block_partition(x,L)

n 	= length(x);
nblocks = floor(n/L);

out = [];
for i=0:nblocks-1
  out = [out x(i*L+1:(i+1)*L,:)];
end

out = reshape(out(:),L^2,nblocks^2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = invblock_partition(x,n,L)

nblocks = floor(n/L);

buf = reshape(x,L,L*nblocks^2);

out = [];
for i=0:nblocks-1
  out = [out;buf(:,i*L*nblocks+1:(i+1)*L*nblocks)];
end


