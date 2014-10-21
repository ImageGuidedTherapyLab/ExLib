import vtk
 
# write vtk points file
def WriteVTKPoints(self,vtkpoints,OutputFileName):
   # loop over points an store in vtk data structure
   # write in meters
   MillimeterMeterConversion = .001;
   scalevtkPoints = vtk.vtkPoints()
   vertices= vtk.vtkCellArray()
   for idpoint in range(vtkpoints.GetNumberOfPoints()):
       point = MillimeterMeterConversion * numpy.array(vtkpoints.GetPoint(idpoint))
       vertices.InsertNextCell( 1 ); vertices.InsertCellPoint( scalevtkPoints.InsertNextPoint(point) )
       #vertices.InsertNextCell( 1 ); vertices.InsertCellPoint( idpoint )

   # set polydata
   polydata = vtk.vtkPolyData()
   polydata.SetPoints(scalevtkPoints )
   polydata.SetVerts( vertices )

   # write to file
   print "WriteVTKPoints: writing",OutputFileName
   polydatawriter = vtk.vtkDataSetWriter()
   polydatawriter.SetFileName(OutputFileName)
   polydatawriter.SetInput(polydata)
   polydatawriter.Update()


DiffusingTipLength  =  100.
DiffusingTipRadius  =  10.

vtkCylinder = vtk.vtkCylinderSource()
vtkCylinder.SetHeight(DiffusingTipLength ); 
vtkCylinder.SetRadius(DiffusingTipRadius );
vtkCylinder.SetCenter(0.0, 0.0, 0.0);
vtkCylinder.SetResolution(16);
vtkCylinder.Update()

# create 3d model
trianglefilter = vtk.vtkDataSetTriangleFilter()
trianglefilter.SetInput(vtkCylinder.GetOutput() )

## cylinderdata = vtkCylinder.GetOutput()
## numpoints = cylinderdata.GetNumberOfPoints()
## 
## # create 3d surface
## # The points to be triangulated are generated randomly in the unit
## # cube located at the origin. The points are then associated with a
## # vtkPolyData.
## math = vtk.vtkMath()
## points = vtk.vtkPoints()
## pointcounter=0
## for iii in range(numpoints ):
##   (xx,yy,zz) =  cylinderdata.GetPoint(iii)
##   for  jjj in range(1,6):
##     scalefactor = jjj/5.
##     points.InsertPoint(pointcounter, scalefactor*xx,scalefactor*yy,scalefactor*zz)
##     pointcounter = pointcounter +1
##     #points.InsertPoint(i, math.Random(0, 1), math.Random(0, 1),
##     #                   math.Random(0, 1))
## 
## profile = vtk.vtkPolyData()
## profile.SetPoints(points)
## 
## # Delaunay3D is used to triangulate the points. The Tolerance is the
## # distance that nearly coincident points are merged
## # together. (Delaunay does better if points are well spaced.) The
## # alpha value is the radius of circumcircles, circumspheres. Any mesh
## # entity whose circumcircle is smaller than this value is output.
## delny = vtk.vtkDelaunay3D()
## delny.SetInput( profile )
## delny.SetTolerance(0.01)
## delny.SetAlpha(0.2)
## delny.BoundingTriangulationOff()
## delny.Update()

# read source landmarks
SourceLMReader = vtk.vtkPolyDataReader()
SourceLMReader.SetFileName("DefaultLandmarks.vtk");
SourceLMReader.Update()
#print SourceLMReader.GetOutput().GetPoints()

# read target landmarks
TargetLMReader = vtk.vtkPolyDataReader()
TargetLMReader.SetFileName("TargetLandmarks.vtk");
TargetLMReader.Update()
#print TargetLMReader.GetOutput().GetPoints()


# create transformation
LandmarkTransform = vtk.vtkLandmarkTransform()
LandmarkTransform.SetSourceLandmarks(SourceLMReader.GetOutput().GetPoints() )
LandmarkTransform.SetTargetLandmarks(TargetLMReader.GetOutput().GetPoints() )
LandmarkTransform.SetModeToRigidBody()
LandmarkTransform.Update()
print LandmarkTransform.GetMatrix()

# apply transform
transformFilter = vtk.vtkTransformFilter()
#transformFilter.SetInput(vtkCylinder.GetOutput() ) 
transformFilter.SetInput(trianglefilter.GetOutput() ) 
transformFilter.SetTransform( LandmarkTransform) 
transformFilter.Update()

# write model
modelWriter = vtk.vtkDataSetWriter()
modelWriter.SetInput(transformFilter.GetOutput())
modelWriter.SetFileName("needle.vtk")
modelWriter.SetFileTypeToBinary()
modelWriter.Update()

# read image/ROI
ImageReader = vtk.vtkDataSetReader()
ImageReader.SetFileName("newimage.vtk")
ImageReader.Update()

# resample needle to image to create mask
vtkResample = vtk.vtkCompositeDataProbeFilter()
vtkResample.SetInput( ImageReader.GetOutput() )
vtkResample.SetSource( transformFilter.GetOutput() ) 
vtkResample.Update()

# write mask
MaskWriter = vtk.vtkDataSetWriter()
MaskWriter.SetInput(vtkResample.GetOutput())
MaskWriter.SetFileName("needlemask.vtk")
MaskWriter.SetFileTypeToBinary()
MaskWriter.Update()

