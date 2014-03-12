figure(1)
set(1,'DoubleBuffer','on');
set(1,'Color',[1 1 1]);

%fid=fopen('tris.bin','rb','ieee-be');
fid=fopen('tris.bin','rb','ieee-le');
tris=fread(fid,inf,'int');
fclose(fid);
tris=reshape(tris,[],3)+1;

%fid=fopen('verts.bin','rb','ieee-be');
fid=fopen('verts.bin','rb','ieee-le');
verts=fread(fid,inf,'float');
fclose(fid);
verts=reshape(verts,[],3);

%disp('Writing surface as OutSurf.vtk');
%write_vtk('OutSurf.vtk',tris,verts);

mys=trisurf(tris,verts(:,1),verts(:,2),verts(:,3),verts(:,1));
  amb=.4; diff=.4; spec=.2;
  set(mys,'FaceColor','interp','EdgeColor','none',...
    'FaceLighting','gour','AmbientStrength',amb,...
    'DiffuseStrength',diff,...
    'SpecularStrength',spec,'EdgeLighting','gour')

axis equal
axis vis3d
view(20,-70);
camlight('headlight');
cameratoolbar
rotate3d on
xlabel('R'); ylabel('A'); zlabel('S');
axis off
shading faceted
colormap(jet)
drawnow

