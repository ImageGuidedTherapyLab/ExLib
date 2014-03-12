function mys=ts(tris,verts)
mys=trisurf(tris,verts(:,1),verts(:,2),verts(:,3));
axis('equal');
rotate3d on
return
