function write_vtk(filename,tris,verts)
%function write_vtk(filename,tris,verts)

% open file
fid=fopen(filename,'w');

% write header
fprintf(fid,'# vtk DataFile Version 2.0\n');
fprintf(fid,'vtk output\n');
fprintf(fid,'ASCII\n');
fprintf(fid,'DATASET POLYDATA\n');

fprintf(fid,'POINTS %d float\n',size(verts,1));
fprintf(fid,'%.03f %.03f %.03f\n',verts');

% write triangles
fprintf(fid,'POLYGONS %d %d\n',size(tris,1),4*size(tris,1));
fprintf(fid,'3 %d %d %d\n',...
  [[tris(:,1)-1] [tris(:,2)-1] [tris(:,3)-1]]');

fclose(fid);
