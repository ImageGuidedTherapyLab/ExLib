M=100;
N=150;
P=50;

%fid=fopen('outvol.bin','rb','ieee-be');
fid=fopen('outvol.bin','rb','ieee-le');
I=fread(fid,[M*N*P],'uint16=>uint16');
fclose(fid);
I=reshape(I,[N M P]);
I=permute(I,[2 1 3]);
mx=max(I,[],3);
k=find(max(mx,[],1));
j=find(max(mx,[],2));    
j=[1 M];
k=[1 N];

figure(2)
set(2,'DoubleBuffer','on');

for ii=1:P
  if length(find(I(:,:,ii)>0))
    w=I(:,:,ii);
    figure(2);
    image(w(min(j):max(j),min(k):max(k)),'CDataMapping','scaled');
    caxis([0 4]);
    axis('image');
    axis('tight');
    colormap([0 0 0; 1 1 1; 0 0 1; 0 1 0; 1 0 0]);
    drawnow
    disp('Paused.  Press a key.');
    pause;
  end
end

