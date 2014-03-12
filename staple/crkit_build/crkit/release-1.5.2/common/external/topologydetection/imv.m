function imv(I,fign)

if nargin>1
  myfig=figure(fign);
else
  myfig=figure;
end

set(myfig,'Color',[1 1 1]);
mn=min(I(:));
mx=max(I(:));
if mx<=mn
  mx=mn+1;
end

image(I,'CDataMapping','scaled');
caxis([mn mx]);
axis('image');
axis('tight');
colormap(gray(256));

return
