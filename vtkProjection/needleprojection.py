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
transformFilter.SetInput(vtkCylinder.GetOutput() ) 
transformFilter.SetTransform( LandmarkTransform) 
transformFilter.Update()

# write model
modelWriter = vtk.vtkDataSetWriter()
modelWriter.SetInput(transformFilter.GetOutput())
modelWriter.SetFileName("needle.vtk")
modelWriter.SetFileTypeToBinary()
modelWriter.Update()

# write model
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
MaskWriter.Update()

