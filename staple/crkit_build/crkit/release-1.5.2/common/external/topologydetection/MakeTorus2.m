M=100;  % rows
N=150;  % cols
P=50;   % planes

I=zeros(M,N,P,'uint16');

I(:,:,round(P/2))=1;

idx=randperm(M*N);
I(idx(1:round(M*N*.5))+(round(P/2)-1)*M*N)=0;

for ii=1:P
  imv(I(:,:,ii),1);
  drawnow
end


fid=fopen('Torus.bin','wb','ieee-be');
for ii=1:P
  fwrite(fid,I(:,:,ii)','uint16');
end
fclose(fid);
