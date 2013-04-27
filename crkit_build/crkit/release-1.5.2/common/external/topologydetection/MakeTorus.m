M=100;  % rows
N=150;  % cols
P=50;   % planes

I=zeros(M,N,P,'uint16');

m=round([M/3 M/3*2]);
n=round([N/3 N/3*2]);
p=round([P/3 P/3*2]);

I(m(1):m(2),n(1):n(2),p(1):p(2))=1;

m=round([M/2-M/12 M/2+M/12]);
n=round([N/2-N/12 N/2+N/12]);
%p=round([P/2-P/12 P/2+P/12]);

I(m(1):m(2),n(1):n(2),:)=0;

for ii=1:P
  imv(I(:,:,ii),1);
  drawnow
end


fid=fopen('Torus.bin','wb','ieee-be');
for ii=1:P
  fwrite(fid,I(:,:,ii)','uint16');
end
fclose(fid);
