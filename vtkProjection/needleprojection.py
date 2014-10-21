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

modelWriter = vtk.vtkDataSetWriter()
modelWriter.SetInput(vtkCylinder.GetOutput())
modelWriter.SetFileName("needle.vtk")
modelWriter.SetFileTypeToBinary()
modelWriter.Update()
